###############################################################################
# Indirect genetic effects (duo model), UNAFFECTED OFFSPRING ONLY:
#   perinatal ~ PRS_child + PRS_mother
###############################################################################

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr)
  library(stringr); library(purrr); library(tibble)
  library(sandwich); library(lmtest)
})

## ---------------------------------------------------------------------------
## CONFIG — adjust paths / column names here
## ---------------------------------------------------------------------------

PRS_DIR        <- "/Users/asgelz01/Downloads/PRS_all/PRS"
PC_DIR         <- "/Users/asgelz01/Downloads/PRS_all"   
META_PATH      <- "/Users/asgelz01/Downloads/SPARK.iWES_v3.2024_08.sample_metadata.tsv"
BASIC_MED_PATH <- "/Users/asgelz01/Downloads/SPARKDataRelease_2026-06-25/basic_medical_screening-2026-06-25.csv"
ROLES_PATH     <- "/Users/asgelz01/Downloads/SPARKDataRelease_2026-06-25/roles-2026-06-25.csv"
INDIV_PATH     <- "/Users/asgelz01/Downloads/SPARKDataRelease_2026-06-25/individuals_registration-2026-06-25.csv"


COL_SPID   <- "spid"
COL_SFID   <- "sfid"
COL_FATHER <- "father"
COL_MOTHER <- "mother"
COL_ASD    <- "asd"
COL_SEX    <- "sex"
COL_BATCH  <- "batch"
COL_TWINS  <- "identical_twins"   # co-twin SPID or TRUE-like flag

N_PCS        <- 10
MIN_EVENTS   <- 5
BETA_SEP_LIM <- 10
TARGET_ANCES <- c("AFR", "AMR", "EUR")

outcomes <- c("birth_etoh_subst", "birth_ivh", "birth_oxygen", "birth_pg_inf",
              "birth_prem", "growth_low_wt", "growth_macroceph", "growth_microceph",
              "med_cond_birth", "med_cond_birth_def", "med_cond_growth")

RESULTS_DIR <- file.path("/Users/asgelz01/Downloads/PRS_all", "ige_duo_regression_unaffected_only")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_PATH <- file.path(RESULTS_DIR, "ige_duo_regression_unaffected_only.tsv")

## ---------------------------------------------------------------------------
## 1. PRS: load all .sscore files, parse trait + ancestry from filename
## ---------------------------------------------------------------------------

sscore_files <- list.files(PRS_DIR, pattern = "\\.sscore$", full.names = TRUE)
if (length(sscore_files) == 0) stop("No .sscore files found in: ", PRS_DIR)

prs_long <- purrr::map_dfr(sscore_files, function(fp) {
  b  <- basename(fp)
  m1 <- stringr::str_match(b, "^([A-Za-z0-9]+)_(AFR|AMR)\\.sscore$")
  m2 <- stringr::str_match(b, "^EUR_([A-Za-z0-9]+)\\.sscore$")
  if (!is.na(m1[1, 1])) { trait <- m1[1, 2]; anc <- m1[1, 3] }
  else if (!is.na(m2[1, 1])) { trait <- m2[1, 2]; anc <- "EUR" }
  else return(NULL)
  if (!anc %in% TARGET_ANCES) return(NULL)
  dt <- data.table::fread(fp, sep = "\t", header = TRUE, check.names = FALSE)
  data.table::setnames(dt, sub("^#", "", names(dt)))
  score_col <- intersect(c("SCORE1_AVG", "SCORE1_SUM", "SCORE"), names(dt))[1]
  if (is.na(score_col)) stop("No score column found in: ", b)
  tibble::tibble(spid  = as.character(dt$IID),
                 trait = trait, ancestry = anc,
                 prs_raw = as.numeric(dt[[score_col]]))
})

prs_traits <- sort(unique(prs_long$trait))
message("PRS traits found: ", paste(prs_traits, collapse = ", "))

# Ancestry assignment per person (from which ancestry-specific file they sit in)
anc_map <- prs_long %>% distinct(spid, ancestry)
dup_anc <- anc_map %>% count(spid) %>% filter(n > 1)
if (nrow(dup_anc) > 0) {
  warning(nrow(dup_anc), " individuals appear in >1 ancestry's PRS files; ",
          "they are DROPPED to keep strata clean.")
  anc_map <- anc_map %>% filter(!spid %in% dup_anc$spid)
  prs_long <- prs_long %>% filter(!spid %in% dup_anc$spid)
}

## ---------------------------------------------------------------------------
## 2. V20 family relationships, offspring sample, sex, ASD, batch
## ---------------------------------------------------------------------------

blank_to_na <- function(x) {
  x <- trimws(as.character(x))
  dplyr::if_else(x %in% c("", "NA", "0", "."), NA_character_, x)
}

split_ids <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", ".")] <- NA_character_
  vals <- x[!is.na(x)]
  if (!length(vals)) return(character(0))
  ids <- trimws(unlist(strsplit(vals, "\\|")))
  unique(ids[nzchar(ids)])
}

parse_asd <- function(x) {
  x <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    x %in% c("true","1","yes","y") ~ 1L,
    x %in% c("false","0","no","n") ~ 0L,
    TRUE ~ NA_integer_
  )
}

roles <- data.table::fread(ROLES_PATH, header = TRUE) %>%
  tibble::as_tibble() %>%
  transmute(
    spid = as.character(subject_sp_id),
    affected_sibs = blank_to_na(affected_sibling_sp_id),
    control_sibs  = blank_to_na(control_sibling_sp_id)
  )

affected_sibling_ids <- split_ids(roles$affected_sibs)
control_sibling_ids  <- split_ids(roles$control_sibs)
child_ids <- unique(c(roles$spid, affected_sibling_ids, control_sibling_ids))


indiv <- data.table::fread(INDIV_PATH, header = TRUE) %>%
  tibble::as_tibble() %>%
  transmute(
    spid = as.character(subject_sp_id),
    sfid = as.character(family_sf_id),
    mother = blank_to_na(biomother_sp_id),
    sex = as.character(sex),
    asd = parse_asd(asd)
  ) %>% distinct(spid, .keep_all = TRUE)

missing_sibling_ids <- setdiff(unique(c(affected_sibling_ids, control_sibling_ids)), indiv$spid)
message("Sibling IDs absent from individuals_registration: ", length(missing_sibling_ids))

meta_genomic <- data.table::fread(META_PATH, sep = "\t", quote = "", header = TRUE) %>%
  tibble::as_tibble() %>%
  transmute(
    spid = as.character(.data[[COL_SPID]]),
    batch = as.character(.data[[COL_BATCH]]),
    twins_raw = trimws(as.character(.data[[COL_TWINS]])),
    twin_partner = dplyr::if_else(stringr::str_detect(twins_raw, "^SP\\d+$"),
                                  twins_raw, NA_character_)
  ) %>% distinct(spid, .keep_all = TRUE)

ped <- indiv %>%
  filter(spid %in% child_ids) %>%
  left_join(meta_genomic, by = "spid")

message("V20 explicit proband/sibling offspring identified: ", nrow(ped))
message("...with own biomother_sp_id available: ", sum(!is.na(ped$mother)))
message("...with known ASD status: ", sum(!is.na(ped$asd)))
message("...unaffected (asd == 0): ", sum(ped$asd == 0L, na.rm = TRUE))

## ---------------------------------------------------------------------------
## 3. Perinatal outcomes (parent-report about the child)
## ---------------------------------------------------------------------------

med <- data.table::fread(BASIC_MED_PATH, header = TRUE) |> tibble::as_tibble()
id_col_med <- intersect(c("subject_sp_id", "spid", "sp_id"), names(med))[1]
if (is.na(id_col_med)) stop("Could not find subject ID column in basic_medical_screening.")

miss_out <- setdiff(outcomes, names(med))
if (length(miss_out) > 0) stop("Outcomes missing from screening file: ",
                               paste(miss_out, collapse = ", "))

clean_v20_binary <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    is.na(x) | x == ""     ~ 0L,
    x == "1"               ~ 1L,
    x == "na_survey_logic" ~ NA_integer_,
    TRUE                   ~ NA_integer_
  )
}

med <- med %>%
  transmute(spid = as.character(.data[[id_col_med]]),
            across(all_of(outcomes), clean_v20_binary)) %>%
  distinct(spid, .keep_all = TRUE)

## ---------------------------------------------------------------------------
## 4. PCs per ancestry 
## ---------------------------------------------------------------------------

pc_files <- list.files(PC_DIR, pattern = "_pcs\\.eigenvec$", full.names = TRUE)
if (length(pc_files) == 0) stop("No *_pcs.eigenvec files found in: ", PC_DIR)

pc_all <- purrr::map_dfr(pc_files, function(fp) {
  fname <- tolower(basename(fp))
  anc <- dplyr::case_when(
    stringr::str_detect(fname, "^afr_") ~ "AFR",
    stringr::str_detect(fname, "^amr_") ~ "AMR",
    stringr::str_detect(fname, "^eur_") ~ "EUR",
    TRUE ~ NA_character_
  )
  if (is.na(anc)) stop("Could not infer ancestry from filename: ", basename(fp))
  dt <- data.table::fread(fp, header = TRUE)
  data.table::setnames(dt, sub("^#", "", names(dt)))
  pc_cols <- paste0("PC", 1:N_PCS)
  if (!all(pc_cols %in% names(dt))) stop("Fewer than ", N_PCS, " PCs in: ", fp)
  dt %>% tibble::as_tibble() %>%
    transmute(spid = as.character(IID), ancestry = anc,
              across(all_of(pc_cols), as.numeric))
})

dup_pc <- pc_all %>% count(spid) %>% filter(n > 1)
if (nrow(dup_pc) > 0)
  stop("Some IIDs appear in multiple ancestry PC files: ",
       paste(head(dup_pc$spid, 10), collapse = ", "))

## ---------------------------------------------------------------------------
## 5. Build duos: UNAFFECTED children only
## ---------------------------------------------------------------------------


children <- ped %>%
  filter(!is.na(mother), !is.na(asd), asd == 0L) %>%
  mutate(tp_key = dplyr::case_when(
    is.na(twin_partner)  ~ NA_character_,
    spid < twin_partner  ~ paste0(spid, "__", twin_partner),
    TRUE                 ~ paste0(twin_partner, "__", spid))) %>%
  arrange(tp_key, spid) %>%
  group_by(tp_key) %>%
  filter(is.na(tp_key) | dplyr::row_number() == 1L) %>%
  ungroup() %>%
  select(spid, sfid, mother, asd, sex, batch) %>%
  inner_join(anc_map, by = "spid") %>%                          # child ancestry
  inner_join(anc_map %>% rename(mother = spid,
                                anc_mother = ancestry), by = "mother") %>%
  filter(ancestry == anc_mother) %>%                            # same stratum
  select(-anc_mother) %>%
  inner_join(med, by = "spid")                                  # has screening record

message("Unaffected mother-child duos (same ancestry, child has screening data): ",
        nrow(children))
print(count(children, ancestry))

## ---------------------------------------------------------------------------
## 6. Residualize PGS on own 10 PCs + batch, per ancestry x role, then scale
## ---------------------------------------------------------------------------

residualize_prs <- function(ids, anc, trait_name) {
  df <- tibble(spid = unique(ids)) %>%
    inner_join(prs_long %>% filter(trait == trait_name, ancestry == anc),
               by = "spid") %>%
    inner_join(pc_all %>% filter(ancestry == anc) %>% select(-ancestry),
               by = "spid") %>%
    left_join(meta_genomic %>% select(spid, batch), by = "spid")
  if (nrow(df) == 0) return(tibble(spid = character(), prs = numeric()))
  pc_terms <- paste(paste0("PC", 1:N_PCS), collapse = " + ")
  use_batch <- !all(is.na(df$batch)) && dplyr::n_distinct(na.omit(df$batch)) > 1
  fml <- if (use_batch)
    as.formula(paste("prs_raw ~", pc_terms, "+ factor(batch)"))
  else
    as.formula(paste("prs_raw ~", pc_terms))
  df   <- df %>% filter(if (use_batch) !is.na(batch) else TRUE)
  fit  <- lm(fml, data = df)
  df$prs <- as.numeric(scale(resid(fit)))
  df %>% select(spid, prs)
}

## ---------------------------------------------------------------------------
## 7. Fit: outcome ~ PRS_child + PRS_mother + sex + child PCs + child batch,
##    cluster SE on mother ID
## ---------------------------------------------------------------------------

fit_duo <- function(dat, outc) {

  # V20 uses na_survey_logic for outcomes that were not administered.
  # Restrict each model to children with an observed value for that outcome.
  dat2 <- dat %>%
    filter(!is.na(.data[[outc]]))

  y <- dat2[[outc]]
  n_case <- sum(y == 1L)
  n_ctrl <- sum(y == 0L)

  base <- tibble(
    outcome = outc,
    n = nrow(dat2),
    n_case = n_case,
    n_control = n_ctrl,
    n_mothers = dplyr::n_distinct(dat2$mother)
  )

  if (nrow(dat2) == 0)
    return(base %>% mutate(status = "no_observed_outcome"))

  if (n_case < MIN_EVENTS || n_ctrl < MIN_EVENTS)
    return(base %>% mutate(status = "skipped_min_events"))

  m <- tryCatch(
    suppressWarnings(
      glm(
        reformulate(
          c("prs_child", "prs_mother", "factor(sex)",
            paste0("PC", 1:N_PCS), "factor(batch)"),
          response = outc
        ),
        family = binomial,
        data = dat2
      )
    ),
    error = function(e) NULL
  )

  if (is.null(m))
    return(base %>% mutate(status = "fit_error"))

  if (!isTRUE(m$converged))
    return(base %>% mutate(status = "not_converged"))

  ct <- tryCatch(
    lmtest::coeftest(
      m,
      vcov = sandwich::vcovCL(
        m,
        cluster = dat2$mother,
        type = "HC1"
      )
    ),
    error = function(e) NULL
  )

  if (is.null(ct))
    return(base %>% mutate(status = "cluster_vcov_error"))
  rows <- purrr::map_dfr(
    c(child = "prs_child", mother = "prs_mother"),
    function(term) {
      if (!term %in% rownames(ct)) return(NULL)
      tibble(beta = ct[term, "Estimate"], se = ct[term, "Std. Error"],
             z = ct[term, "z value"], p = ct[term, "Pr(>|z|)"])
    }, .id = "coef")
  bind_cols(base[rep(1, nrow(rows)), ], rows) %>%
    mutate(status = dplyr::case_when(
             !is.finite(se) | se <= 0 | se > 10 ~ "unstable_se",
             abs(beta) > BETA_SEP_LIM ~ "separation",
             TRUE ~ "ok"),
           or = exp(beta),
           or_lo = exp(beta - 1.96 * se),
           or_hi = exp(beta + 1.96 * se))
}

results <- purrr::map_dfr(TARGET_ANCES, function(anc) {
  kids <- children %>% filter(ancestry == anc)
  if (nrow(kids) == 0) return(NULL)
  purrr::map_dfr(prs_traits, function(tr) {
    prs_c <- residualize_prs(kids$spid,   anc, tr) %>% rename(prs_child  = prs)
    prs_m <- residualize_prs(kids$mother, anc, tr) %>%
      rename(mother = spid, prs_mother = prs)
    child_pcs <- pc_all %>%
      filter(ancestry == anc) %>%
      select(spid, all_of(paste0("PC", 1:N_PCS)))

    dat <- kids %>%
      inner_join(prs_c, by = "spid") %>%
      inner_join(prs_m, by = "mother") %>%
      inner_join(child_pcs, by = "spid") %>%
      filter(
        !is.na(sex),
        !is.na(batch),
        !if_any(all_of(paste0("PC", 1:N_PCS)), ~ is.na(.x))
      )
    stopifnot(!anyDuplicated(dat$spid))
    if (nrow(dat) == 0) return(NULL)
    purrr::map_dfr(outcomes, ~ fit_duo(dat, .x)) %>%
      mutate(ancestry = anc, trait = tr, .before = 1)
  })
})

## ---------------------------------------------------------------------------
## 8. Stouffer meta across ancestries 
## ---------------------------------------------------------------------------

meta_res <- results %>%
  filter(status == "ok") %>%
  mutate(n_eff = 4 / (1 / n_case + 1 / n_control),
         z_signed = beta / se) %>%
  group_by(trait, outcome, coef) %>%
  filter(dplyr::n_distinct(ancestry) == 3) %>%
  summarise(
    ancestry   = "META",
    k          = dplyr::n(),
    n          = sum(n),
    n_case     = sum(n_case),
    n_control  = sum(n_control),
    n_mothers  = sum(n_mothers),
    z          = sum(sqrt(n_eff) * z_signed) / sqrt(sum(n_eff)),
    p          = 2 * pnorm(-abs(z)),
    beta       = sum(n_eff * beta) / sum(n_eff),   # descriptive pooled beta
    se         = NA_real_, or = exp(beta), or_lo = NA_real_, or_hi = NA_real_,
    status     = "ok",
    .groups = "drop"
  )

out <- bind_rows(results, meta_res) %>%
  group_by(coef, ancestry, trait) %>%
  mutate(q_fdr = dplyr::if_else(status == "ok",
                                p.adjust(replace(p, status != "ok", NA),
                                         method = "BH"), NA_real_)) %>%
  ungroup() %>%
  arrange(coef, ancestry, trait, outcome)

data.table::fwrite(out, OUT_PATH, sep = "\t")
message("Wrote: ", OUT_PATH, "  (", nrow(out), " rows)")

## ---------------------------------------------------------------------------
## 9. PAPER-READY TABLES 
## ---------------------------------------------------------------------------

# Labels matching the main pipeline's trait_map / outcome_map
TRAIT_LABELS <- c(AS = "Asthma", OB = "Obesity", SCZ = "Schizophrenia",
                  MDD = "Major depression", PPH = "Postpartum hemorrhage",
                  PE = "Pre-eclampsia", GD = "Gestational diabetes")

OUTCOME_LABELS <- c(
  birth_etoh_subst   = "Fetal alcohol syndrome or in-utero drug/alcohol exposure",
  birth_ivh          = "Intraventricular hemorrhage",
  birth_oxygen       = "Oxygen deprivation at birth requiring NICU stay",
  birth_pg_inf       = "Serious prenatal infection",
  birth_prem         = "Premature birth",
  growth_low_wt      = "Difficulty gaining weight",
  growth_macroceph   = "Macrocephaly",
  growth_microceph   = "Microcephaly",
  med_cond_birth     = "Birth or pregnancy complications",
  med_cond_birth_def = "Birth defects",
  med_cond_growth    = "Growth abnormalities")

TERM_LABELS <- c(child  = "Child PGS (direct)",
                 mother = "Maternal PGS (indirect)")

num2  <- function(x) formatC(x, digits = 2, format = "f")
p_fmt <- function(p) ifelse(p < 1e-4, formatC(p, format = "e", digits = 1),
                            formatC(p, digits = 3, format = "f"))
q_fmt <- p_fmt

build_paper_table <- function(res_sub) {
  res_sub %>%
    filter(coef %in% names(TERM_LABELS)) %>%
    mutate(
      trait_full   = dplyr::coalesce(TRAIT_LABELS[trait], trait),
      outcome_full = dplyr::coalesce(OUTCOME_LABELS[outcome], outcome),
      Term         = factor(TERM_LABELS[coef], levels = unname(TERM_LABELS)),
      sig_bh     = !is.na(q_fdr) & q_fdr < 0.05,
      `Effect (OR)` = dplyr::if_else(
        is.na(or), "NA",
        dplyr::if_else(is.na(or_lo), num2(or),
                       paste0(num2(or), " [", num2(or_lo), "\u2013",
                              num2(or_hi), "]"))),
      Display = paste0(
        "\u03b2 ", ifelse(is.na(beta), "NA", num2(beta)),
        ifelse(is.na(se), "", paste0(" (", num2(se), ")")), "; ",
        "OR ", ifelse(is.na(or), "NA", num2(or)),
        ifelse(is.na(or_lo), "",
               paste0(" [", num2(or_lo), "\u2013", num2(or_hi), "]")),
        "; p = ", ifelse(is.na(p), "NA", p_fmt(p)),
        "; BH q = ", ifelse(is.na(q_fdr), "NA", q_fmt(q_fdr)),
        ifelse(sig_bh, " *", ""))
    ) %>%
    arrange(trait_full, outcome_full, Term) %>%
    select(
      `Trait`           = trait_full,
      `Outcome`         = outcome_full,
      `Term`,
      `Effect (OR)`,
      `Beta`            = beta,
      `Standard error`  = se,
      `Odds ratio`      = or,
      `CI lower`        = or_lo,
      `CI upper`        = or_hi,
      `p-value`         = p,
      `BH q-value`      = q_fdr,
      `BH significant`  = sig_bh,
      `Individuals (n)` = n,
      `Mothers (n)`     = n_mothers,
      `Cases`           = n_case,
      `Controls`        = n_control,
      Display
    )
}


for (anc in TARGET_ANCES) {
  sub <- out %>% filter(ancestry == anc, status == "ok")
  if (nrow(sub) == 0) next

  fp <- file.path(RESULTS_DIR, paste0("ige_duo_unaffected_paper_table_", anc, ".tsv"))
  data.table::fwrite(build_paper_table(sub), fp, sep = "\t")
  message("Wrote paper table: ", fp)
}

meta_sub <- out %>% filter(ancestry == "META", status == "ok")

if (nrow(meta_sub) > 0) {
  meta_paper <- meta_sub %>%
    filter(coef %in% names(TERM_LABELS)) %>%
    mutate(
      trait_full   = dplyr::coalesce(TRAIT_LABELS[trait], trait),
      outcome_full = dplyr::coalesce(OUTCOME_LABELS[outcome], outcome),
      Term         = factor(TERM_LABELS[coef], levels = unname(TERM_LABELS)),
      sig_bh       = !is.na(q_fdr) & q_fdr < 0.05,
      `Effect (OR)` = ifelse(is.na(or), "NA", num2(or)),
      Display = paste0(
        "OR ", ifelse(is.na(or), "NA", num2(or)),
        "; Meta Z = ", ifelse(is.na(z), "NA", num2(z)),
        "; p = ", ifelse(is.na(p), "NA", p_fmt(p)),
        "; BH q = ", ifelse(is.na(q_fdr), "NA", q_fmt(q_fdr)),
        ifelse(sig_bh, " *", "")
      )
    ) %>%
    arrange(trait_full, outcome_full, Term) %>%
    select(
      `Trait`           = trait_full,
      `Outcome`         = outcome_full,
      `Term`,
      `Effect (OR)`,
      `Odds ratio`      = or,
      `Meta Z`          = z,
      `p-value`         = p,
      `BH q-value`      = q_fdr,
      `BH significant`  = sig_bh,
      `Individuals (n)` = n,
      `Mothers (n)`     = n_mothers,
      `Cases`           = n_case,
      `Controls`        = n_control,
      Display
    )

  meta_fp <- file.path(RESULTS_DIR, "ige_duo_unaffected_paper_table_META.tsv")
  data.table::fwrite(meta_paper, meta_fp, sep = "\t")
  message("Wrote META paper table: ", meta_fp)
}

## Sanity summary
message("\nUnaffected duo counts by ancestry:")
print(results %>% distinct(ancestry, trait, n) %>%
        group_by(ancestry) %>% summarise(max_n = max(n)))
message("\nFDR-significant (q < .05), meta:")
print(out %>% filter(ancestry == "META", q_fdr < 0.05) %>%
        select(coef, trait, outcome, or, p, q_fdr))
