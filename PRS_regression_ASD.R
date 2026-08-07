# =====================================================
# PRS ~ Perinatal Outcomes
# Ancestry-specific regressions for all traits
# Stouffer meta-analysis for selected cross-ancestry traits
# =====================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(stringr)
})

# -----------------------------
# Paths
# -----------------------------
PRS_DIR        <- "/Users/asgelz01/Downloads/PRS_all/PRS_ASD"
PC_DIR         <- "/Users/asgelz01/Downloads/PRS_all"
BASIC_MED_PATH <- "/Users/asgelz01/Downloads/SPARKDataRelease_2025-07-14/basic_medical_screening-2025-07-14.csv"
META_PATH      <- "/Users/asgelz01/Downloads/SPARK.iWES_v3.2024_08.sample_metadata.tsv"

OUT_DIR <- file.path(PRS_DIR, "per_ancestry_outputs")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Outcomes
# -----------------------------
outcomes_bin <- c(
  "birth_etoh_subst","birth_ivh","birth_oxygen","birth_pg_inf","birth_prem",
  "growth_low_wt","growth_macroceph","growth_microceph",
  "med_cond_birth","med_cond_birth_def","med_cond_growth"
)

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
  med_cond_growth    = "Growth abnormalities"
)

# -----------------------------
# Load phenotype data
# -----------------------------
basic_med <- fread(BASIC_MED_PATH) %>%
  select(subject_sp_id, all_of(outcomes_bin)) %>%
  rename(IID = subject_sp_id) %>%
  mutate(
    across(all_of(outcomes_bin), ~ tidyr::replace_na(., 0L)),
    across(all_of(outcomes_bin), as.integer)
  )

# -----------------------------
# Load covariates
# -----------------------------
meta_cov <- fread(META_PATH) %>%
  transmute(
    IID   = spid,
    sex   = as.factor(sex),
    batch = as.factor(batch)
  ) %>%
  distinct(IID, .keep_all = TRUE)

# -----------------------------
# Load ancestry-specific PCs
# -----------------------------
pc_files <- list.files(PC_DIR, pattern = "_pcs\\.eigenvec$", full.names = TRUE)
if (length(pc_files) == 0) stop("No *_pcs.eigenvec files found in: ", PC_DIR)

infer_ancestry <- function(path) {
  fname <- tolower(basename(path))
  case_when(
    str_detect(fname, "^afr_") ~ "AFR",
    str_detect(fname, "^amr_") ~ "AMR",
    str_detect(fname, "^eur_") ~ "EUR",
    TRUE ~ NA_character_
  )
}

read_pcs <- function(f) {
  anc <- infer_ancestry(f)
  if (is.na(anc)) stop("Could not infer ancestry from filename: ", basename(f))
  
  dt <- fread(f)
  pc_cols <- grep("^PC\\d+$", names(dt), value = TRUE)
  if (length(pc_cols) == 0) stop("No PC columns in: ", f)
  
  dt %>%
    select(IID, all_of(pc_cols)) %>%
    mutate(
      across(all_of(pc_cols), as.numeric),
      ancestry_code = anc
    )
}

pcs <- map_dfr(pc_files, read_pcs)

dup_iids <- pcs %>%
  count(IID) %>%
  filter(n > 1) %>%
  pull(IID)

if (length(dup_iids) > 0) {
  stop(
    "Some IIDs appear in multiple ancestry PC files: ",
    paste(head(dup_iids, 10), collapse = ", ")
  )
}

pc_names <- grep("^PC\\d+$", names(pcs), value = TRUE)

# -----------------------------
# Load PRS .sscore files
# -----------------------------
read_sscore <- function(f) {
  df <- fread(f)
  score_col <- if ("SCORESUM" %in% names(df)) "SCORESUM" else names(df)[ncol(df)]
  
  df %>%
    transmute(
      IID = IID,
      PRS = as.numeric(.data[[score_col]]),
      PRS_name = tools::file_path_sans_ext(basename(f))
    )
}

sscore_files <- list.files(PRS_DIR, pattern = "\\.sscore$", full.names = TRUE)
sscore_files <- sscore_files[!grepl("^META_", basename(sscore_files))]

if (length(sscore_files) == 0) {
  stop("No non-META .sscore files found in: ", PRS_DIR)
}

prs_wide <- map_dfr(sscore_files, read_sscore) %>%
  pivot_wider(names_from = PRS_name, values_from = PRS)

# -----------------------------
# Merge data
# -----------------------------
dat0 <- basic_med %>%
  inner_join(meta_cov, by = "IID") %>%
  inner_join(pcs, by = "IID") %>%
  inner_join(prs_wide, by = "IID")

prs_names <- setdiff(
  names(dat0),
  c("IID", outcomes_bin, "sex", "batch", "ancestry_code", pc_names)
)

# -----------------------------
# Fit ancestry-specific logistic models
# -----------------------------
fit_one_ancestry <- function(df_anc, anc_label) {
  
  df_anc <- df_anc %>%
    mutate(across(all_of(prs_names), ~ as.numeric(scale(.x))))
  
  covar_terms <- c("sex", "batch", pc_names)
  covar_rhs <- paste(covar_terms, collapse = " + ")
  
  fit_logistic <- function(outcome, prs) {
    
    fml <- as.formula(paste0(outcome, " ~ ", prs, " + ", covar_rhs))
    
    df2 <- df_anc %>%
      filter(
        !is.na(.data[[outcome]]),
        !if_any(all_of(c(pc_names, prs)), ~ is.na(.x))
      )
    
    if (nrow(df2) < 50) return(NULL)
    
    n_cases <- sum(df2[[outcome]] == 1, na.rm = TRUE)
    n_controls <- sum(df2[[outcome]] == 0, na.rm = TRUE)
    
    if (n_cases == 0 || n_controls == 0) return(NULL)
    
    n_eff <- 4 / ((1 / n_cases) + (1 / n_controls))
    
    fit <- tryCatch(
      glm(fml, data = df2, family = binomial()),
      error = function(e) NULL
    )
    
    if (is.null(fit)) return(NULL)
    
    sm <- summary(fit)$coefficients
    if (!(prs %in% rownames(sm))) return(NULL)
    
    beta <- sm[prs, "Estimate"]
    se   <- sm[prs, "Std. Error"]
    p    <- sm[prs, "Pr(>|z|)"]
    
    tibble(
      ancestry = anc_label,
      outcome = outcome,
      PRS = prs,
      model = "logistic",
      estimate = beta,
      std_error = se,
      p_value = p,
      OR = exp(beta),
      CI_low = exp(beta - 1.96 * se),
      CI_high = exp(beta + 1.96 * se),
      n = nobs(fit),
      n_cases = n_cases,
      n_controls = n_controls,
      n_eff = n_eff
    )
  }
  
  map_dfr(outcomes_bin, function(outcome) {
    map_dfr(prs_names, function(prs) {
      fit_logistic(outcome, prs)
    })
  })
}

# -----------------------------
# Trait extraction and labels
# -----------------------------
trait_map <- c(
  AS  = "Asthma",
  OB  = "Obesity",
  SCZ = "Schizophrenia",
  MDD = "Major depression",
  PPH = "Postpartum hemorrhage",
  PE  = "Pre-eclampsia",
  GD  = "Gestational diabetes"
)

traits_keep <- c(
  "Schizophrenia",
  "Obesity",
  "Asthma",
  "Major depression"
)

anc_map_full <- c(
  EUR = "European",
  AFR = "African",
  AMR = "Admixed American"
)

# -----------------------------
# Run per-ancestry regressions
# -----------------------------
for (code in c("EUR", "AFR", "AMR")) {
  
  anc_label <- anc_map_full[[code]]
  df_anc <- dat0 %>% filter(ancestry_code == code)
  
  if (nrow(df_anc) == 0) {
    message("Skipping ", anc_label, " because no rows were found.")
    next
  }
  
  message("Running ancestry: ", anc_label, " n=", nrow(df_anc))
  
  res_anc <- fit_one_ancestry(df_anc, anc_label)
  
  if (is.null(res_anc) || nrow(res_anc) == 0) {
    message("No results for ", anc_label)
    next
  }
  
  res_anc <- res_anc %>%
    mutate(
      trait_code = case_when(
        str_starts(PRS, "EUR_")  ~ str_split_fixed(PRS, "_", 3)[, 2],
        str_starts(PRS, "META_") ~ str_split_fixed(PRS, "_", 3)[, 2],
        TRUE                     ~ str_split_fixed(PRS, "_", 3)[, 1]
      ),
      trait = recode(trait_code, !!!trait_map, .default = trait_code)
    ) %>%
    group_by(trait) %>%
    mutate(
      q_BH = p.adjust(p_value, method = "BH")
    ) %>%
    ungroup() %>%
    relocate(trait, .after = PRS) %>%
    select(-trait_code)
  
  out_path <- file.path(
    OUT_DIR,
    paste0("PRS_outcome_", gsub(" ", "_", anc_label), ".tsv")
  )
  
  fwrite(res_anc, out_path, sep = "\t")
  message("Wrote: ", out_path)
}


# =========================================
# Paper-ready tables 
# =========================================

# --- Helpers ---
p_fmt <- function(x) ifelse(is.na(x), "",
                            ifelse(x < 1e-6, formatC(x, format="e", digits=2),
                                   ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))))
q_fmt <- p_fmt
num2  <- function(x) ifelse(is.na(x), "", sprintf("%.2f", x))

read_ancestry_tsv <- function(anc_label) {
  f <- file.path(OUT_DIR, paste0("PRS_outcome_", gsub(" ", "_", anc_label), ".tsv"))
  if (!file.exists(f)) {
    message("No file found for ", anc_label, ": ", f)
    return(NULL)
  }
  readr::read_tsv(f, show_col_types = FALSE)
}

prep_results_tbl <- function(df, anc_label) {
  if (is.null(df) || !nrow(df)) return(NULL)
  
  # Ensure BH exists (per-trait within-ancestry)
  if (!"q_BH" %in% names(df)) {
    df <- df %>%
      dplyr::group_by(trait) %>%
      dplyr::mutate(q_BH = p.adjust(p_value, "BH")) %>%
      dplyr::ungroup()
  }

  df <- df %>%
    dplyr::mutate(
      ancestry      = anc_label,
      trait         = dplyr::recode(trait, !!!trait_map, .default = trait),
      outcome_label = dplyr::recode(outcome, !!!outcome_map, .default = outcome),

      OR       = suppressWarnings(as.numeric(OR)),
      CI_lower = suppressWarnings(as.numeric(CI_low)),
      CI_upper = suppressWarnings(as.numeric(CI_high))
    )
  
  df %>%
    dplyr::select(trait, outcome_label,
                  estimate, std_error, p_value, q_BH,
                  OR, CI_lower, CI_upper,
                  n
    )
}

# ==============================================
# Paper table
# ==============================================
write_paper_table <- function(df, anc_label, filename) {
  if (is.null(df) || !nrow(df)) return(invisible(NULL))
  
  out <- df %>%
    dplyr::group_by(trait) %>%
    dplyr::arrange(dplyr::desc(abs(estimate)), .by_group = TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      `Effect (OR)` = paste0(num2(OR), " [", num2(CI_lower), "–", num2(CI_upper), "]"),
      Display = paste0(
        "β ", num2(estimate), " (", num2(std_error), "); ",
        "OR ", num2(OR), " [", num2(CI_lower), "–", num2(CI_upper), "]; ",
        "p = ", p_fmt(p_value), "; BH q = ", q_fmt(q_BH),
        ifelse(!is.na(q_BH) && q_BH < 0.05, " *", "")
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(`Significant (BH q<0.05)` = q_BH < 0.05) %>%
    dplyr::select(
      `Trait`                 = trait,
      `Outcome`               = outcome_label,
      `Effect (OR)`,
      `Beta`                  = estimate,       # log-odds coefficient
      `Standard error`        = std_error,
      `CI lower`              = CI_lower,
      `CI upper`              = CI_upper,
      `p-value`               = p_value,
      `BH q-value`            = q_BH,
      `Significant (BH q<0.05)`,
      `Sample size (n)`       = n,
      `Display`
    )
  
  readr::write_csv(out, filename)
  message("Wrote: ", filename)
  return(out)
}


eur_raw <- read_ancestry_tsv("European")
afr_raw <- read_ancestry_tsv("African")
amr_raw <- read_ancestry_tsv("Admixed American")

eur_tbl <- prep_results_tbl(eur_raw, "European")
afr_tbl <- prep_results_tbl(afr_raw, "African")
amr_tbl <- prep_results_tbl(amr_raw, "Admixed American")

write_paper_table(eur_tbl, "European",         file.path(OUT_DIR, "prs_perinatal_European.csv"))
write_paper_table(afr_tbl, "African",          file.path(OUT_DIR, "prs_perinatal_African.csv"))
write_paper_table(amr_tbl, "Admixed American", file.path(OUT_DIR, "prs_perinatal_Admixed_American.csv"))


# =====================================================
# Effective-sample-size weighted Stouffer meta-analysis
# =====================================================

files <- c(
  European = file.path(OUT_DIR, "PRS_outcome_European.tsv"),
  African = file.path(OUT_DIR, "PRS_outcome_African.tsv"),
  `Admixed American` = file.path(OUT_DIR, "PRS_outcome_Admixed_American.tsv")
)

res_all <- imap_dfr(files, function(path, anc_name) {
  if (!file.exists(path)) stop("Missing file: ", path)
  
  read_tsv(path, show_col_types = FALSE) %>%
    mutate(ancestry = anc_name)
})

required_cols <- c(
  "ancestry", "trait", "outcome", "estimate", "std_error",
  "p_value", "n", "n_cases", "n_controls", "n_eff"
)

missing_cols <- setdiff(required_cols, names(res_all))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

res_all_clean <- res_all %>%
  filter(trait %in% traits_keep) %>%
  mutate(
    outcome_label = recode(outcome, !!!outcome_map, .default = outcome),
    estimate = as.numeric(estimate),
    std_error = as.numeric(std_error),
    p_value = as.numeric(p_value),
    n = as.numeric(n),
    n_cases = as.numeric(n_cases),
    n_controls = as.numeric(n_controls),
    n_eff = as.numeric(n_eff),
    
    z_ancestry = estimate / std_error,
    weight_stouffer = sqrt(n_eff)
  ) %>%
  filter(
    !is.na(estimate),    is.finite(estimate),    abs(estimate) <= 10,
    !is.na(std_error),   is.finite(std_error),   std_error > 0,  std_error <= 10,
    !is.na(z_ancestry),  is.finite(z_ancestry),
    !is.na(weight_stouffer), is.finite(weight_stouffer), weight_stouffer > 0,
    !is.na(n_eff),       is.finite(n_eff),        n_eff > 0,
    n_cases > 0,
    n_controls > 0
  )

# Keep only trait × outcome combinations with valid estimates in all three ancestry strata 
res_all_clean <- res_all_clean %>%
  group_by(trait, outcome_label) %>%
  filter(n_distinct(ancestry) == 3) %>%
  ungroup()


write_csv(
  res_all_clean,
  file.path(OUT_DIR, "stouffer_neff_weighted_meta_input_selected_traits.csv")
)

meta_stouffer <- res_all_clean %>%
  group_by(trait, outcome_label) %>%
  summarise(
    z_meta = sum(weight_stouffer * z_ancestry, na.rm = TRUE) /
      sqrt(sum(weight_stouffer^2, na.rm = TRUE)),
    
    p_value_meta = 2 * pnorm(abs(z_meta), lower.tail = FALSE),
    
    n_total = sum(n, na.rm = TRUE),
    n_cases_total = sum(n_cases, na.rm = TRUE),
    n_controls_total = sum(n_controls, na.rm = TRUE),
    n_eff_total = sum(n_eff, na.rm = TRUE),
    
    k_ancestries = n_distinct(ancestry),
    
    ancestries_included = paste(
      sort(unique(ancestry)),
      collapse = "; "
    ),
    
    beta_neff_weighted = sum(n_eff * estimate, na.rm = TRUE) /
      sum(n_eff, na.rm = TRUE),
    
    OR_neff_weighted = exp(beta_neff_weighted),
    
    .groups = "drop"
  ) %>%
  group_by(trait) %>%
  mutate(
    q_BH_meta = p.adjust(p_value_meta, method = "BH")
  ) %>%
  ungroup() %>%
  arrange(trait, q_BH_meta, p_value_meta)

write_csv(
  meta_stouffer,
  file.path(OUT_DIR, "stouffer_neff_weighted_meta_selected_traits.csv")
)

print(meta_stouffer)

# =========================================
# Paper-ready table for Stouffer meta-analysis
# =========================================

meta_paper <- meta_stouffer %>%
  mutate(
    `Effect (OR)` = paste0(
      num2(OR_neff_weighted),
      " (descriptive; no CI†)"
    ),
    Display = paste0(
      "β = ", num2(beta_neff_weighted), "; ",
      "OR ", num2(OR_neff_weighted), " (descriptive); ",
      "Z = ", num2(z_meta), "; ",
      "p = ", p_fmt(p_value_meta), "; ",
      "BH q = ", q_fmt(q_BH_meta),
      ifelse(!is.na(q_BH_meta) & q_BH_meta < 0.05, " *", "")
    ),
    
    `Significant (BH q<0.05)` = q_BH_meta < 0.05
  ) %>%
  select(
    `Trait` = trait,
    `Outcome` = outcome_label,
    
    `Meta Z` = z_meta,
    `Meta p-value` = p_value_meta,
    `Meta BH q-value` = q_BH_meta,
    `Significant (BH q<0.05)`,
    
    `Beta, N_eff-weighted` = beta_neff_weighted,
    `OR, N_eff-weighted` = OR_neff_weighted,
    `Effect (OR)`,
    
    `N total` = n_total,
    `Cases total` = n_cases_total,
    `Controls total` = n_controls_total,
    `Effective sample size total` = n_eff_total,
    
    `Number of ancestries` = k_ancestries,
    `Ancestries included` = ancestries_included,
    
    `Display`
  )

write_csv(
  meta_paper,
  file.path(
    OUT_DIR,
    "prs_perinatal_stouffer_meta_paper_ready.csv"
  )
)

print(meta_paper)
