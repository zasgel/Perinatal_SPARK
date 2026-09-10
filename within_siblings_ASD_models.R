# =========================================================
# Within-siblings PRS analysis — ASD status models — SPARK V20
#
#   MODEL = "asd_cov"  y ~ PRS_mean + PRS_within + asd + sex + batch + (1|FID)
#   MODEL = "asd_int"  ... + PRS_within:asd
#
# =========================================================

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(tibble)
  library(lme4); library(purrr); library(stringr); library(readr)
})

# -----------------------------
# Paths
# -----------------------------
SSCORE_DIR     <- "/Users/asgelz01/Downloads/PRS_all/PRS"
META_PATH      <- "/Users/asgelz01/Downloads/SPARK.iWES_v3.2024_08.sample_metadata.tsv"
BASIC_MED_PATH <- "/Users/asgelz01/Downloads/SPARKDataRelease_2026-06-25/basic_medical_screening-2026-06-25.csv"
ROLES_PATH     <- "/Users/asgelz01/Downloads/SPARKDataRelease_2026-06-25/roles-2026-06-25.csv"
INDIV_PATH     <- "/Users/asgelz01/Downloads/SPARKDataRelease_2026-06-25/individuals_registration-2026-06-25.csv"

# -----------------------------
# Toggles
# -----------------------------
MODEL      <- "asd_int"   # "asd_cov" | "asd_int"
CHILD_ONLY <- TRUE     
ANCESTRIES <- c("EUR", "AFR", "AMR")
MIN_EVENTS <- 5           
DIGITS_MZ  <- 6          
BH_SCOPE <- "within_trait"

OUT_DIR <- file.path("/Users/asgelz01/Downloads/PRS_all",
                     paste0("within_sibs_", MODEL,
                            if (CHILD_ONLY) "_childonly" else "_allmembers"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Labels
# -----------------------------
trait_map <- c(AS = "Asthma", OB = "Obesity", SCZ = "Schizophrenia",
               MDD = "Major depression", PPH = "Postpartum hemorrhage",
               PE = "Pre-eclampsia", GD = "Gestational diabetes")

outcome_map <- c(
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

term_map <- c(
  "PRS_mean"       = "Between-family (PRS_mean)",
  "PRS_within"     = "Within-family (PRS_within)",
  "asd"            = "ASD status",
  "PRS_within:asd" = "Within-family PRS × ASD interaction")

term_order <- c("PRS_mean", "PRS_within", "asd", "PRS_within:asd")

outcomes     <- names(outcome_map)
anc_map_full <- c(EUR = "European", AFR = "African", AMR = "Admixed American")

TERMS_KEEP <- if (MODEL == "asd_cov") {
  c("PRS_mean", "PRS_within", "asd")
} else {
  c("PRS_mean", "PRS_within", "asd", "PRS_within:asd")
}

# -----------------------------
# Load data
# -----------------------------
blank_na <- function(x) {
  x <- trimws(as.character(x))
  dplyr::if_else(x %in% c("", "NA", "0", "."), NA_character_, x)
}

meta_raw <- fread(META_PATH, sep = "\t") %>%
  as_tibble() %>%
  mutate(spid = as.character(spid))

# -----------------------------
# MZ twin identification from legacy iWES genomic metadata
# -----------------------------
parse_twin_tokens <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "0", ".", "False", "FALSE", "false",
             "True", "TRUE", "true", "1")] <- NA_character_
  strsplit(x, "[,;/|[:space:]]+")
}

build_mz_groups <- function(ids, twin_col) {
  toks <- parse_twin_tokens(twin_col)
  edges <- map2_dfr(ids, toks, function(i, tt) {
    tt <- tt[!is.na(tt) & nzchar(tt) & tt != i]
    if (!length(tt)) return(NULL)
    tibble(a = i, b = tt)
  })
  if (!nrow(edges)) return(tibble(IID = character(0), mz_group = character(0)))

  nodes <- unique(c(edges$a, edges$b))
  parent <- setNames(nodes, nodes)
  find_root <- function(x) {
    while (parent[[x]] != x) {
      parent[[x]] <<- parent[[parent[[x]]]]
      x <- parent[[x]]
    }
    x
  }
  for (k in seq_len(nrow(edges))) {
    ra <- find_root(edges$a[k]); rb <- find_root(edges$b[k])
    if (ra != rb) parent[[rb]] <- ra
  }
  tibble(IID = nodes,
         mz_group = vapply(nodes, find_root, character(1), USE.NAMES = FALSE))
}

mz_groups <- build_mz_groups(meta_raw$spid, meta_raw$identical_twins)
frac_known <- if (nrow(mz_groups)) mean(mz_groups$IID %in% meta_raw$spid) else 0

if (nrow(mz_groups) > 0 && frac_known > 0.5) {
  MZ_MODE <- "spid"
  message("MZ mode: SPID-linked. ", nrow(mz_groups), " individuals in ",
          n_distinct(mz_groups$mz_group), " MZ sets.")
} else {
  MZ_MODE <- "prs_fallback"
  message("MZ mode: FALLBACK (PRS similarity). Inspect identical_twins before trusting fallback.")
}

# -----------------------------
# V20 explicit proband + sibling construction
# -----------------------------
parse_asd <- function(x) {
  x <- tolower(trimws(as.character(x)))
  case_when(
    x %in% c("true","1","yes","y") ~ 1L,
    x %in% c("false","0","no","n") ~ 0L,
    TRUE ~ NA_integer_
  )
}

split_ids <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", ".")] <- NA_character_
  vals <- x[!is.na(x)]
  if (!length(vals)) return(character(0))
  ids <- trimws(unlist(strsplit(vals, "\\|")))
  unique(ids[nzchar(ids)])
}

roles <- fread(ROLES_PATH) %>%
  as_tibble() %>%
  transmute(
    IID = as.character(subject_sp_id),
    affected_sibs = blank_na(affected_sibling_sp_id),
    control_sibs  = blank_na(control_sibling_sp_id)
  )

affected_sibling_ids <- split_ids(roles$affected_sibs)
control_sibling_ids  <- split_ids(roles$control_sibs)
child_ids <- unique(c(roles$IID, affected_sibling_ids, control_sibling_ids))

indiv <- fread(INDIV_PATH) %>%
  as_tibble() %>%
  transmute(
    IID = as.character(subject_sp_id),
    FID = as.character(family_sf_id),
    sex = factor(sex),
    asd = parse_asd(asd)
  ) %>%
  distinct(IID, .keep_all = TRUE)

missing_sibling_ids <- setdiff(
  unique(c(affected_sibling_ids, control_sibling_ids)), indiv$IID
)
message("Sibling IDs absent from individuals_registration: ", length(missing_sibling_ids))

genomic_cov <- meta_raw %>%
  transmute(
    IID = as.character(spid),
    batch = factor(batch),
    mz_flag = !(identical_twins %in% c(NA, "", "False", "FALSE", "false", "0", "NA"))
  ) %>%
  left_join(mz_groups, by = "IID") %>%
  distinct(IID, .keep_all = TRUE)

meta <- indiv %>%
  filter(IID %in% child_ids, !is.na(FID), !is.na(sex), !is.na(asd)) %>%
  inner_join(genomic_cov, by = "IID") %>%
  select(IID, FID, sex, batch, mz_flag, mz_group, asd) %>%
  distinct(IID, .keep_all = TRUE)

if (!isTRUE(CHILD_ONLY)) {
  stop("CHILD_ONLY must remain TRUE in the V20 within-sibling ASD models.")
}

message("V20 proband/sibling sample with known ASD status: ", nrow(meta),
        " individuals in ", n_distinct(meta$FID), " families | ASD cases: ",
        sum(meta$asd == 1L), " | unaffected: ", sum(meta$asd == 0L))

# -----------------------------
# V20 phenotype coding
# -----------------------------
clean_v20_binary <- function(x) {
  x <- trimws(as.character(x))
  unexpected <- unique(x[
    !is.na(x) & x != "" & !x %in% c("1", "na_survey_logic")
  ])
  if (length(unexpected) > 0) {
    warning("Unexpected V20 outcome value(s): ",
            paste(head(unexpected, 10), collapse = ", "))
  }
  case_when(
    is.na(x) | x == ""     ~ 0L,
    x == "1"               ~ 1L,
    x == "na_survey_logic" ~ NA_integer_,
    TRUE                   ~ NA_integer_
  )
}

basic_med <- fread(BASIC_MED_PATH, colClasses = "character") %>%
  as_tibble() %>%
  select(subject_sp_id, all_of(outcomes)) %>%
  rename(IID = subject_sp_id) %>%
  mutate(
    IID = as.character(IID),
    across(all_of(outcomes), clean_v20_binary)
  ) %>%
  distinct(IID, .keep_all = TRUE)

# -----------------------------
# Formatting helpers
# -----------------------------
num2  <- function(x, d = 2, na = "") ifelse(is.na(x), na, formatC(x, format="f", digits=d))
p_fmt <- function(p, na = "") ifelse(is.na(p), na,
           ifelse(p < 1e-4, formatC(p, format="e", digits=2),
                            formatC(p, format="f", digits=3)))
q_fmt <- p_fmt

# -----------------------------
# Helpers
# -----------------------------
read_sscore_simple <- function(path) {
  df <- fread(path, sep = "\t", header = TRUE, check.names = FALSE)
  setnames(df, sub("^#", "", names(df)))
  if (!("IID" %in% names(df))) stop("IID missing in ", path)
  sc <- if ("SCORESUM" %in% names(df)) "SCORESUM" else
        if ("SCORE1_SUM" %in% names(df)) "SCORE1_SUM" else names(df)[ncol(df)]
  tibble(IID = as.character(df$IID), PRS_raw = as.numeric(df[[sc]]))
}


keep_one_per_mz_group <- function(df) {
  if (!nrow(df)) return(df)
  df %>%
    group_by(FID, mz_group) %>%
    filter(is.na(mz_group) | IID == min(IID)) %>%
    ungroup()
}

keep_one_per_mz_pair_prs <- function(df, digits = DIGITS_MZ) {
  if (!nrow(df)) return(df)
  df$..z <- as.numeric(scale(df$PRS_raw))
  out <- bind_rows(lapply(split(df, df$FID, drop = TRUE), function(d) {
    idx <- which(d$mz_flag %in% TRUE)
    if (length(idx) < 2) return(d)
    grp <- round(d$..z[idx], digits)
    drop_ids <- character(0)
    for (g in unique(grp)) {
      ids <- d$IID[idx][grp == g]
      if (length(ids) > 1) drop_ids <- c(drop_ids, sort(ids)[-1])
    }
    d[!(d$IID %in% drop_ids), , drop = FALSE]
  }))
  out$..z <- NULL
  out
}

dedup_mz <- function(df) {
  if (MZ_MODE == "spid") keep_one_per_mz_group(df) else keep_one_per_mz_pair_prs(df)
}


prep_within_family_data <- function(prs_tbl, meta_anc, basic_med) {
  d <- inner_join(meta_anc, basic_med, by = "IID") %>%
    inner_join(prs_tbl, by = "IID") %>%
    filter(!is.na(FID), !is.na(sex), !is.na(batch)) %>%
    add_count(FID, name = "fam_n") %>% filter(fam_n >= 2) %>% select(-fam_n)

  n_before <- nrow(d)
  d <- dedup_mz(d)
  n_mz_dropped <- n_before - nrow(d)
  if (n_mz_dropped > 0)
    message("    MZ dedup removed ", n_mz_dropped, " individual(s)")

  d %>%
    add_count(FID, name = "fam_n") %>% filter(fam_n >= 2) %>% select(-fam_n) %>%
    mutate(PRS_z = as.numeric(scale(PRS_raw))) %>%
    group_by(FID) %>%
    mutate(PRS_mean = mean(PRS_z, na.rm = TRUE), PRS_within = PRS_z - PRS_mean) %>%
    ungroup()
}

get_safe_vcov <- function(fit) {
  V <- tryCatch(as.matrix(vcov(fit, use.hessian = FALSE)), error = function(e) NULL)
  if (!is.null(V) && all(is.finite(diag(V))) && all(diag(V) >= 0)) return(V)
  V <- tryCatch(as.matrix(vcov(fit)), error = function(e) NULL)
  if (!is.null(V) && all(is.finite(diag(V))) && all(diag(V) >= 0)) return(V)
  NULL
}

build_formula <- function(yvar) {
  rhs <- if (MODEL == "asd_cov") {
    "PRS_mean + PRS_within + asd"
  } else {
    "PRS_mean + PRS_within + asd + PRS_within:asd"
  }
  as.formula(paste0(yvar, " ~ ", rhs, " + sex + batch + (1|FID)"))
}

# -----------------------------
# Model fit
# -----------------------------
fit_one_condition <- function(yvar, d) {
  empty <- function(note, dd = NULL, ev = NA_integer_, nev = NA_integer_,
                    ndisc = NA_integer_, evdisc = NA_integer_,
                    nasd = NA_integer_) tibble(
    condition = yvar, term = NA_character_, estimate = NA_real_, std.error = NA_real_,
    z = NA_real_, p = NA_real_, OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
    n_ind = if (is.null(dd)) 0L else nrow(dd),
    n_fam = if (is.null(dd)) 0L else n_distinct(dd$FID),
    n_disc_fam = ndisc, events_disc = evdisc, n_asd_disc_fam = nasd,
    events = ev, nonevents = nev, note = note)

  if (is.null(d) || !nrow(d)) return(empty("No data"))

  dd <- d %>% filter(!is.na(.data[[yvar]]))
  if (nrow(dd) < 4) return(empty("Too few rows", dd))
  dd <- dd %>% add_count(FID, name = "fam_n") %>% filter(fam_n >= 2) %>% select(-fam_n)
  if (nrow(dd) < 4) return(empty("No families with >=2 members", dd))

  # families where the OUTCOME varies -- these identify PRS_within
  fam <- dd %>% group_by(FID) %>%
    summarise(n1 = sum(.data[[yvar]] == 1L), n0 = sum(.data[[yvar]] == 0L),
              .groups = "drop") %>%
    mutate(disc = n1 > 0 & n0 > 0)
  n_disc  <- sum(fam$disc)
  ev_disc <- sum(fam$n1[fam$disc])

  # families where ASD STATUS varies -- these identify `asd` and the interactions
  n_asd_disc <- dd %>% group_by(FID) %>%
    summarise(v = n_distinct(asd), .groups = "drop") %>% filter(v > 1) %>% nrow()

  events <- sum(dd[[yvar]] == 1, na.rm = TRUE)
  nonev  <- sum(dd[[yvar]] == 0, na.rm = TRUE)

  if (n_disc == 0)
    return(empty("No within-family variation in outcome", dd, events, nonev,
                 n_disc, ev_disc, n_asd_disc))

  if (events < MIN_EVENTS || nonev < MIN_EVENTS)
    return(empty(paste0("Too few events or non-events (< ", MIN_EVENTS, ")"),
                 dd, events, nonev, n_disc, ev_disc, n_asd_disc))
  if (ev_disc < MIN_EVENTS)
    return(empty(paste0("Too few events in discordant families (< ", MIN_EVENTS, ")"),
                 dd, events, nonev, n_disc, ev_disc, n_asd_disc))
  if (n_distinct(dd$asd) < 2)
    return(empty("No ASD variation", dd, events, nonev, n_disc, ev_disc, n_asd_disc))
  if (MODEL == "asd_int" && n_asd_disc == 0)
    return(empty("No ASD-discordant families", dd, events, nonev,
                 n_disc, ev_disc, n_asd_disc))

  fit <- tryCatch(
    glmer(build_formula(yvar), data = dd, family = binomial(),
          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5),
            calc.derivs = TRUE, check.conv.grad = "ignore",
            check.conv.singular = "ignore", check.conv.hess = "ignore"),
          nAGQ = 0),
    error = function(e) e)

  if (inherits(fit, "error"))
    return(empty(paste("Model error:", conditionMessage(fit)), dd, events, nonev,
                 n_disc, ev_disc, n_asd_disc))

  fe <- lme4::fixef(fit)

  stability_terms <- if (MODEL == "asd_cov") {
    intersect(c("PRS_mean", "PRS_within"), names(fe))
  } else {
    intersect(c("PRS_mean", "PRS_within", "PRS_within:asd"), names(fe))
  }

  if (length(stability_terms) == 0)
    return(empty("PRS term(s) missing from fitted model", dd, events, nonev,
                 n_disc, ev_disc, n_asd_disc))

  if (any(abs(fe[stability_terms]) > 10, na.rm = TRUE))
    return(empty("Possible separation/instability in PRS term: |beta| > 10",
                 dd, events, nonev, n_disc, ev_disc, n_asd_disc))

  V <- get_safe_vcov(fit)
  if (!is.null(V)) {
    se <- sqrt(pmax(0, diag(V)))
    names(se) <- names(fe)

    prs_se <- se[stability_terms]
    if (any(!is.finite(prs_se)) || any(prs_se <= 0) || any(prs_se > 10)) {
      return(empty("Unstable PRS term: non-finite, non-positive, or SE > 10",
                   dd, events, nonev, n_disc, ev_disc, n_asd_disc))
    }

    z_val <- ifelse(se > 0, fe / se, NA_real_)
    p_val <- ifelse(is.finite(z_val), 2 * pnorm(abs(z_val), lower.tail = FALSE), NA_real_)
    OR <- exp(fe); CI_low <- exp(fe - 1.96*se); CI_high <- exp(fe + 1.96*se)
    note_msg <- NA_character_
  } else {
    se <- z_val <- p_val <- OR <- CI_low <- CI_high <- rep(NA_real_, length(fe))
    note_msg <- "Non-PD vcov"
  }

  tibble(condition = yvar, term = names(fe), estimate = as.numeric(fe),
         std.error = as.numeric(se), z = as.numeric(z_val), p = as.numeric(p_val),
         OR = as.numeric(OR), CI_low = as.numeric(CI_low), CI_high = as.numeric(CI_high),
         n_ind = nrow(dd), n_fam = n_distinct(dd$FID),
         n_disc_fam = n_disc, events_disc = ev_disc, n_asd_disc_fam = n_asd_disc,
         events = events, nonevents = nonev, note = note_msg) %>%
    filter(term %in% TERMS_KEEP)
}

# -----------------------------
# BH correction 
# -----------------------------
apply_bh <- function(res) {
  if (BH_SCOPE == "within_trait") {
    res %>% group_by(trait_full, term) %>%
      mutate(q_BH = p.adjust(p, method = "BH")) %>% ungroup()
  } else if (BH_SCOPE == "across_traits") {
    res %>% group_by(term) %>%
      mutate(q_BH = p.adjust(p, method = "BH")) %>% ungroup()
  } else stop("Unknown BH_SCOPE: ", BH_SCOPE)
}

# -----------------------------
# Paper-ready table builder (ALL terms in one table)
# -----------------------------
make_paper_table <- function(res) {
  res %>%
    filter(!is.na(estimate)) %>%
    mutate(
      term_label = recode(term, !!!term_map, .default = term),
      term_rank  = match(term, term_order),
      `Effect (OR)` = paste0(num2(OR), " [", num2(CI_low), "-", num2(CI_high), "]"),
      Display = paste0("OR ", num2(OR), " [", num2(CI_low), "-", num2(CI_high),
                       "]; p = ", p_fmt(p), "; BH q = ", q_fmt(q_BH),
                       ifelse(!is.na(q_BH) & q_BH < 0.05, " *", "")),
      `BH significant` = !is.na(q_BH) & q_BH < 0.05) %>%
    arrange(ancestry, trait_full, outcome_label, term_rank) %>%
    select(Ancestry = ancestry, Trait = trait_full, Outcome = outcome_label,
           Term = term_label, `Effect (OR)`,
           Beta = estimate, `Standard error` = std.error,
           `Odds ratio` = OR, `CI lower` = CI_low, `CI upper` = CI_high,
           `p-value` = p, `BH q-value` = q_BH, `BH significant`,
           `Individuals (n)` = n_ind, `Families (n)` = n_fam,
           `Outcome-discordant families` = n_disc_fam,
           `Events in discordant families` = events_disc,
           `ASD-discordant families` = n_asd_disc_fam,
           Events = events, `Non-events` = nonevents, Display)
}

# -----------------------------
# PRS file discovery
# -----------------------------
files_all <- list.files(SSCORE_DIR, pattern = "\\.sscore$", full.names = TRUE)

label_file <- function(fp) {
  b <- tools::file_path_sans_ext(basename(fp))
  parts <- str_split_fixed(b, "_", 2); p1 <- parts[1]; p2 <- parts[2]
  anc_codes <- c("EUR","AFR","AMR"); trait_codes <- names(trait_map)
  if (p1 %in% anc_codes && p2 %in% trait_codes) return(tibble(ancestry=p1, trait=p2, file=fp))
  if (p1 %in% trait_codes && p2 %in% anc_codes) return(tibble(ancestry=p2, trait=p1, file=fp))
  NULL
}
labels <- map_dfr(files_all, label_file)
if (nrow(labels) == 0) stop("No PRS files matched expected naming patterns.")

# -----------------------------
# Run per ancestry
# -----------------------------
all_anc_res <- list()

for (code in ANCESTRIES) {
  anc_label  <- anc_map_full[[code]]
  anc_label_file <- gsub(" ", "_", anc_label)
  labels_anc <- labels %>% filter(ancestry == code)
  if (nrow(labels_anc) == 0) { message("No PRS files for ", code); next }

  message("\n=== ", anc_label, " | model=", MODEL,
          " | child_only=", CHILD_ONLY, " | mz_mode=", MZ_MODE,
          " | n_files=", nrow(labels_anc), " ===")

  res_anc_list <- list()
  for (tr in unique(labels_anc$trait)) {
    prs_tbl <- map_dfr(labels_anc$file[labels_anc$trait == tr], read_sscore_simple) %>%
      distinct(IID, .keep_all = TRUE)
    d_tr <- prep_within_family_data(prs_tbl, meta, basic_med)
    message(sprintf("  trait=%-4s n=%d in %d families (%d ASD-discordant)",
                    tr, nrow(d_tr), n_distinct(d_tr$FID),
                    d_tr %>% group_by(FID) %>%
                      summarise(v = n_distinct(asd), .groups="drop") %>%
                      filter(v > 1) %>% nrow()))
    res_anc_list[[tr]] <- map_dfr(outcomes, ~ fit_one_condition(.x, d_tr)) %>%
      mutate(trait = tr, ancestry = anc_label)
  }

  res_anc <- bind_rows(res_anc_list) %>%
    mutate(trait_full = recode(trait, !!!trait_map, .default = trait),
           outcome_label = recode(condition, !!!outcome_map, .default = condition)) %>%
    apply_bh()

  # raw results
  fwrite(res_anc, file.path(OUT_DIR,
    paste0("within_sibs_", MODEL, "_", anc_label_file, ".tsv")), sep = "\t")

  # paper-ready table for THIS ancestry, all terms together
  paper_anc <- make_paper_table(res_anc)
  write_csv(paper_anc, file.path(OUT_DIR,
    paste0("paper_ready_", MODEL, "_", anc_label_file, ".csv")))

  message("  wrote ", anc_label, ": ", nrow(res_anc), " rows raw, ",
          nrow(paper_anc), " rows paper-ready")
  all_anc_res[[anc_label]] <- res_anc
}

res_all <- bind_rows(all_anc_res)
fwrite(res_all, file.path(OUT_DIR,
  paste0("within_sibs_", MODEL, "_all_ancestries.tsv")), sep = "\t")

# combined paper-ready across all ancestries (q_BH kept from per-ancestry runs)
write_csv(make_paper_table(res_all),
          file.path(OUT_DIR, paste0("paper_ready_", MODEL, "_all_ancestries.csv")))

# -----------------------------
# Summary
# -----------------------------
cat("\n===== MODEL:", MODEL, "| CHILD_ONLY:", CHILD_ONLY,
    "| MZ_MODE:", MZ_MODE, "| BH_SCOPE:", BH_SCOPE, "=====\n")
cat("\nFitted models by term:\n")
print(as.data.frame(res_all %>% filter(!is.na(estimate)) %>%
                      count(ancestry, term, name = "n_models")))

cat("\nDropped models, by reason:\n")
print(as.data.frame(res_all %>% filter(is.na(estimate)) %>% count(note, sort = TRUE)))

if (MODEL == "asd_int") {
  it <- "PRS_within:asd"
  cat("\nInteraction: ", it, "\n", sep = "")
  tb <- res_all %>% filter(term == it, !is.na(estimate)) %>%
    select(ancestry, trait_full, outcome_label, n_fam, n_disc_fam,
           n_asd_disc_fam, estimate, std.error, p, q_BH) %>% arrange(p)
  if (nrow(tb)) print(as.data.frame(tb)) else cat("  none fitted\n")
} else {
  cat("\nASD covariate:\n")
  tb <- res_all %>% filter(term == "asd", !is.na(estimate)) %>%
    select(ancestry, trait_full, outcome_label, n_fam, n_asd_disc_fam,
           estimate, std.error, p, q_BH) %>% arrange(p)
  if (nrow(tb)) print(as.data.frame(head(tb, 20))) else cat("  none fitted\n")
}

message("\nWrote results to: ", OUT_DIR)
