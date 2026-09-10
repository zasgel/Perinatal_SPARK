###############################################################################
# Indirect genetic effects (duo model): perinatal ~ PRS_child * ASD_child + PRS_mother, interaction
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
PC_DIR         <- "/Users/asgelz01/Downloads/PRS_all"   # {afr,amr,eur}_pcs.eigenvec live here
META_PATH      <- "/Users/asgelz01/Downloads/SPARK.iWES_v3.2024_08.sample_metadata.tsv"
BASIC_MED_PATH <- "/Users/asgelz01/Downloads/SPARKDataRelease_2025-07-14/basic_medical_screening-2025-07-14.csv"

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

RESULTS_DIR <- file.path("/Users/asgelz01/Downloads/PRS_all", "ige_duo_regression_ASD_interaction")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_PATH <- file.path(RESULTS_DIR, "ige_duo_regression_ASD_interaction.tsv")

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

anc_map <- prs_long %>% distinct(spid, ancestry)
dup_anc <- anc_map %>% count(spid) %>% filter(n > 1)
if (nrow(dup_anc) > 0) {
  warning(nrow(dup_anc), " individuals appear in >1 ancestry's PRS files; ",
          "they are DROPPED to keep strata clean.")
  anc_map <- anc_map %>% filter(!spid %in% dup_anc$spid)
  prs_long <- prs_long %>% filter(!spid %in% dup_anc$spid)
}

## ---------------------------------------------------------------------------
## 2. Metadata: pedigree, ASD status, sex, batch
## ---------------------------------------------------------------------------

meta <- data.table::fread(META_PATH, sep = "\t", quote = "", header = TRUE) |>
  tibble::as_tibble()

need_cols <- c(COL_SPID, COL_SFID, COL_FATHER, COL_MOTHER, COL_ASD, COL_SEX,
               COL_TWINS)
if (!is.na(COL_BATCH)) need_cols <- c(need_cols, COL_BATCH)
missing_cols <- setdiff(need_cols, names(meta))
if (length(missing_cols) > 0) {
  stop("Metadata is missing column(s): ", paste(missing_cols, collapse = ", "),
       "\nAvailable columns:\n  ", paste(names(meta), collapse = "\n  "))
}

blank_to_na <- function(x) {
  x <- trimws(as.character(x))
  dplyr::if_else(x %in% c("", "NA", "0", "."), NA_character_, x)
}

ped <- meta %>%
  transmute(
    spid   = as.character(.data[[COL_SPID]]),
    sfid   = as.character(.data[[COL_SFID]]),
    father = blank_to_na(.data[[COL_FATHER]]),
    mother = blank_to_na(.data[[COL_MOTHER]]),
    asd    = dplyr::case_when(
      tolower(trimws(as.character(.data[[COL_ASD]]))) %in% c("true","1","yes","y") ~ 1L,
      tolower(trimws(as.character(.data[[COL_ASD]]))) %in% c("false","0","no","n") ~ 0L,
      TRUE ~ NA_integer_
    ),
    sex    = as.character(.data[[COL_SEX]]),
    batch  = if (!is.na(COL_BATCH)) as.character(.data[[COL_BATCH]]) else NA_character_,
    twins_raw = trimws(as.character(.data[[COL_TWINS]]))
  ) %>%
  mutate(twin_partner = dplyr::if_else(stringr::str_detect(twins_raw, "^SP\\d+$"),
                                       twins_raw, NA_character_))

parent_ids <- unique(na.omit(c(ped$father, ped$mother)))
ped <- ped %>%
  mutate(has_parent = !is.na(father) | !is.na(mother),
         is_parent  = spid %in% parent_ids,
         is_child   = has_parent & !is_parent)

message("Children identified in pedigree: ", sum(ped$is_child))

## ---------------------------------------------------------------------------
## 3. Perinatal outcomes (parent-report about the child)
## ---------------------------------------------------------------------------

med <- data.table::fread(BASIC_MED_PATH, header = TRUE) |> tibble::as_tibble()
id_col_med <- intersect(c("subject_sp_id", "spid", "sp_id"), names(med))[1]
if (is.na(id_col_med)) stop("Could not find subject ID column in basic_medical_screening.")

miss_out <- setdiff(outcomes, names(med))
if (length(miss_out) > 0) stop("Outcomes missing from screening file: ",
                               paste(miss_out, collapse = ", "))

med <- med %>%
  transmute(spid = as.character(.data[[id_col_med]]),
            across(all_of(outcomes),
                   ~ dplyr::if_else(is.na(as.numeric(.x)), 0, as.numeric(.x)))) %>%
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
## 5. Build duos
## ---------------------------------------------------------------------------

children <- ped %>%
  filter(is_child, !is.na(mother), !is.na(asd)) %>%
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

message("Mother-child duos (same ancestry, child has screening data): ",
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
    left_join(ped %>% select(spid, batch), by = "spid")
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
## 7. Fit: outcome ~ PRS_child * ASD_child + PRS_mother + sex + child PCs + child batch,
##    cluster SE on mother ID
## ---------------------------------------------------------------------------

fit_duo <- function(dat, outc) {
  y <- dat[[outc]]
  n_case <- sum(y == 1); n_ctrl <- sum(y == 0)
  base <- tibble(outcome = outc, n = nrow(dat),
                 n_case = n_case, n_control = n_ctrl,
                 n_mothers = dplyr::n_distinct(dat$mother),
                 n_asd = sum(dat$asd == 1, na.rm = TRUE),
                 n_unaffected = sum(dat$asd == 0, na.rm = TRUE))
  if (n_case < MIN_EVENTS || n_ctrl < MIN_EVENTS)
    return(base %>% mutate(status = "skipped_min_events"))

  m <- tryCatch(
    suppressWarnings(
      glm(
        reformulate(
          c("prs_child * asd", "prs_mother", "factor(sex)",
            paste0("PC", 1:N_PCS), "factor(batch)"),
          response = outc
        ),
        family = binomial,
        data = dat
      )
    ),
    error = function(e) NULL
  )
  if (is.null(m))
    return(base %>% mutate(status = "fit_error"))
  if (!m$converged)
    return(base %>% mutate(status = "not_converged"))

  ct <- lmtest::coeftest(m, vcov = sandwich::vcovCL(m, cluster = dat$mother,
                                                    type = "HC1"))
  rows <- purrr::map_dfr(
    c(
      child_unaffected = "prs_child",
      mother = "prs_mother",
      child_x_asd = "prs_child:asd"
    ),
    function(term) {
      if (!term %in% rownames(ct)) return(NULL)
      tibble(beta = ct[term, "Estimate"], se = ct[term, "Std. Error"],
             z = ct[term, "z value"], p = ct[term, "Pr(>|z|)"])
    }, .id = "coef")
  bind_cols(base[rep(1, nrow(rows)), ], rows) %>%
    mutate(status = dplyr::if_else(abs(beta) > BETA_SEP_LIM,
                                   "separation", "ok"),
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
        !is.na(asd),
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
## 8. Ancestry-specific multiple-testing correction and outputs
## ---------------------------------------------------------------------------

out <- results %>%
  group_by(coef, ancestry, trait) %>%
  mutate(
    q_fdr = dplyr::if_else(
      status == "ok",
      p.adjust(replace(p, status != "ok", NA), method = "BH"),
      NA_real_
    )
  ) %>%
  ungroup() %>%
  arrange(coef, ancestry, trait, outcome)

data.table::fwrite(out, OUT_PATH, sep = "\t")
message("Wrote: ", OUT_PATH, "  (", nrow(out), " rows)")

## ---------------------------------------------------------------------------
## 9. PAPER-READY TABLES — ancestry-specific only
## ---------------------------------------------------------------------------

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

TERM_LABELS <- c(
  child_unaffected = "Child PGS effect (unaffected, ASD=0)",
  mother           = "Maternal PGS (conditional)",
  child_x_asd      = "Child PGS x ASD interaction"
)

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
      sig_bh       = !is.na(q_fdr) & q_fdr < 0.05,
      `Effect (OR)` = dplyr::if_else(
        is.na(or), "NA",
        dplyr::if_else(
          is.na(or_lo), num2(or),
          paste0(num2(or), " [", num2(or_lo), "–", num2(or_hi), "]")
        )
      ),
      Display = paste0(
        "β ", ifelse(is.na(beta), "NA", num2(beta)),
        ifelse(is.na(se), "", paste0(" (", num2(se), ")")), "; ",
        "OR ", ifelse(is.na(or), "NA", num2(or)),
        ifelse(is.na(or_lo), "",
               paste0(" [", num2(or_lo), "–", num2(or_hi), "]")),
        "; p = ", ifelse(is.na(p), "NA", p_fmt(p)),
        "; BH q = ", ifelse(is.na(q_fdr), "NA", q_fmt(q_fdr)),
        ifelse(sig_bh, " *", "")
      )
    ) %>%
    arrange(trait_full, outcome_full, Term) %>%
    select(
      `Trait`                    = trait_full,
      `Outcome`                  = outcome_full,
      `Term`,
      `Effect (OR)`,
      `Beta`                     = beta,
      `Standard error`           = se,
      `Odds ratio`               = or,
      `CI lower`                 = or_lo,
      `CI upper`                 = or_hi,
      `p-value`                  = p,
      `BH q-value`               = q_fdr,
      `BH significant`           = sig_bh,
      `Individuals (n)`          = n,
      `Mothers (n)`              = n_mothers,
      `ASD children (n)`         = n_asd,
      `Unaffected children (n)`  = n_unaffected,
      `Cases`                    = n_case,
      `Controls`                 = n_control,
      Display
    )
}

for (anc in TARGET_ANCES) {
  sub <- out %>% filter(ancestry == anc, status == "ok")
  if (nrow(sub) == 0) next

  fp <- file.path(
    RESULTS_DIR,
    paste0("ige_duo_regression_ASD_interaction_paper_table_", anc, ".tsv")
  )
  data.table::fwrite(build_paper_table(sub), fp, sep = "\t")
  message("Wrote paper table: ", fp)
}

message("\nDuo counts by ancestry:")
print(
  results %>%
    distinct(ancestry, trait, n) %>%
    group_by(ancestry) %>%
    summarise(max_n = max(n), .groups = "drop")
)

message("\nFDR-significant child/maternal/interaction effects (q < .05), by ancestry:")
print(
  out %>%
    filter(status == "ok", q_fdr < 0.05) %>%
    select(ancestry, coef, trait, outcome, or, p, q_fdr)
)

message("\nNo cross-ancestry meta-analysis performed for the ASD-interaction duo model.")
