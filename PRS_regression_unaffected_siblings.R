# =====================================================
# PRS ~ Perinatal Outcomes  —  POPULATION-LEVEL, UNAFFECTED SIBLINGS
# =====================================================

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(purrr)
  library(readr); library(stringr); library(sandwich); library(lmtest)
})

# -----------------------------
# Paths
# -----------------------------
PRS_DIR        <- "/Users/asgelz01/Downloads/PRS_all/PRS"  
PC_DIR         <- "/Users/asgelz01/Downloads/PRS_all"
BASIC_MED_PATH <- "/Users/asgelz01/Downloads/SPARKDataRelease_2026-06-25/basic_medical_screening-2026-06-25.csv"
ROLES_PATH     <- "/Users/asgelz01/Downloads/SPARKDataRelease_2026-06-25/roles-2026-06-25.csv"
INDIV_PATH     <- "/Users/asgelz01/Downloads/SPARKDataRelease_2026-06-25/individuals_registration-2026-06-25.csv"
IWES_META_PATH <- "/Users/asgelz01/Downloads/SPARK.iWES_v3.2024_08.sample_metadata.tsv"

# -----------------------------
# Toggles
# -----------------------------
CLUSTER_SE   <- TRUE   
CLUSTER_TYPE <- "HC1"   

OUT_DIR <- file.path(PC_DIR, "population_clustered_unaffected_siblings")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Outcomes
# -----------------------------
outcomes_bin <- c(
  "birth_etoh_subst","birth_ivh","birth_oxygen","birth_pg_inf","birth_prem",
  "growth_low_wt","growth_macroceph","growth_microceph",
  "med_cond_birth","med_cond_birth_def","med_cond_growth")

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

trait_map <- c(AS = "Asthma", OB = "Obesity", SCZ = "Schizophrenia",
               MDD = "Major depression", PPH = "Postpartum hemorrhage",
               PE = "Pre-eclampsia", GD = "Gestational diabetes")

traits_keep  <- c("Schizophrenia", "Obesity", "Asthma", "Major depression")
anc_map_full <- c(EUR = "European", AFR = "African", AMR = "Admixed American")

# -----------------------------
# Phenotypes — V20 coding
# -----------------------------
clean_v20_binary <- function(x) {
  x <- trimws(as.character(x))
  case_when(
    is.na(x) | x == ""     ~ 0L,
    x == "1"               ~ 1L,
    x == "na_survey_logic" ~ NA_integer_,
    TRUE                   ~ NA_integer_
  )
}

basic_med <- fread(BASIC_MED_PATH, colClasses = "character") %>%
  as_tibble() %>%
  select(subject_sp_id, all_of(outcomes_bin)) %>%
  rename(IID = subject_sp_id) %>%
  mutate(IID = as.character(IID),
         across(all_of(outcomes_bin), clean_v20_binary)) %>%
  distinct(IID, .keep_all = TRUE)

# -----------------------------
# V20 explicit proband/sibling sample; restrict to unaffected (ASD=0)
# -----------------------------
blank_na <- function(x) {
  x <- trimws(as.character(x))
  if_else(x %in% c("", "NA", "0", "."), NA_character_, x)
}
parse_asd <- function(x) {
  x <- tolower(trimws(as.character(x)))
  case_when(x %in% c("true","1","yes","y") ~ 1L,
            x %in% c("false","0","no","n") ~ 0L,
            TRUE ~ NA_integer_)
}
split_ids <- function(x) {
  x <- trimws(as.character(x)); x[x %in% c("", "NA", ".")] <- NA_character_
  vals <- x[!is.na(x)]
  if (!length(vals)) return(character(0))
  ids <- trimws(unlist(strsplit(vals, "\\|")))
  unique(ids[nzchar(ids)])
}

roles <- fread(ROLES_PATH) %>% as_tibble() %>%
  transmute(IID=as.character(subject_sp_id),
            biomother=blank_na(biomother_sp_id),
            biofather=blank_na(biofather_sp_id),
            affected_sibs=blank_na(affected_sibling_sp_id),
            control_sibs=blank_na(control_sibling_sp_id))

affected_sibling_ids <- split_ids(roles$affected_sibs)
control_sibling_ids  <- split_ids(roles$control_sibs)
child_ids <- unique(c(roles$IID, affected_sibling_ids, control_sibling_ids))

indiv <- fread(INDIV_PATH) %>% as_tibble() %>%
  transmute(IID=as.character(subject_sp_id),
            FID=as.character(family_sf_id),
            sex=factor(sex), asd=parse_asd(asd)) %>%
  distinct(IID,.keep_all=TRUE)

missing_sibling_ids <- setdiff(unique(c(affected_sibling_ids,control_sibling_ids)), indiv$IID)
message("Sibling IDs absent from individuals_registration: ", length(missing_sibling_ids))

meta_cov <- indiv %>%
  filter(IID %in% child_ids, !is.na(FID), !is.na(asd), asd == 0L) %>%
  distinct(IID,.keep_all=TRUE)

batch_cov <- fread(IWES_META_PATH, sep="\t") %>% as_tibble() %>%
  transmute(IID=as.character(spid), batch=factor(batch)) %>%
  distinct(IID,.keep_all=TRUE)

message("Unaffected V20 sample: ", nrow(meta_cov), " individuals in ",
        n_distinct(meta_cov$FID), " families.")

# -----------------------------
# Ancestry-specific PCs
# -----------------------------
pc_files <- list.files(PC_DIR, pattern = "_pcs\\.eigenvec$", full.names = TRUE)
if (length(pc_files) == 0) stop("No *_pcs.eigenvec files found in: ", PC_DIR)

infer_ancestry <- function(path) {
  fname <- tolower(basename(path))
  case_when(str_detect(fname, "^afr_") ~ "AFR",
            str_detect(fname, "^amr_") ~ "AMR",
            str_detect(fname, "^eur_") ~ "EUR",
            TRUE ~ NA_character_)
}

read_pcs <- function(f) {
  anc <- infer_ancestry(f)
  if (is.na(anc)) stop("Could not infer ancestry from filename: ", basename(f))
  dt <- fread(f)
  pc_cols <- grep("^PC\\d+$", names(dt), value = TRUE)
  if (length(pc_cols) == 0) stop("No PC columns in: ", f)
  dt %>% select(IID, all_of(pc_cols)) %>%
    mutate(IID = as.character(IID),
           across(all_of(pc_cols), as.numeric), ancestry_code = anc)
}

pcs <- map_dfr(pc_files, read_pcs)
dup_iids <- pcs %>% count(IID) %>% filter(n > 1) %>% pull(IID)
if (length(dup_iids) > 0)
  stop("IIDs in multiple ancestry PC files: ", paste(head(dup_iids, 10), collapse = ", "))
pc_names <- grep("^PC\\d+$", names(pcs), value = TRUE)

# -----------------------------
# PRS
# -----------------------------
read_sscore <- function(f) {
  df <- fread(f)
  setnames(df, sub("^#", "", names(df)))
  score_col <- if ("SCORESUM" %in% names(df)) "SCORESUM" else
               if ("SCORE1_SUM" %in% names(df)) "SCORE1_SUM" else names(df)[ncol(df)]
  df %>% transmute(IID = as.character(IID),
                   PRS = as.numeric(.data[[score_col]]),
                   PRS_name = tools::file_path_sans_ext(basename(f)))
}

sscore_files <- list.files(PRS_DIR, pattern = "\\.sscore$", full.names = TRUE)
if (length(sscore_files) == 0) stop("No .sscore files found in: ", PRS_DIR)

prs_wide <- map_dfr(sscore_files, read_sscore) %>%
  pivot_wider(names_from = PRS_name, values_from = PRS)

# -----------------------------
# Merge
# -----------------------------
dat0 <- basic_med %>%
  inner_join(meta_cov, by = "IID") %>%
  inner_join(batch_cov, by = "IID") %>%
  inner_join(pcs, by = "IID") %>%
  inner_join(prs_wide, by = "IID")

prs_names <- setdiff(names(dat0),
  c("IID","FID", outcomes_bin, "sex","batch","asd", "ancestry_code", pc_names))

message("\nAnalytic sample: ", nrow(dat0), " unaffected individuals in ",
        n_distinct(dat0$FID), " families")
print(dat0 %>% count(ancestry_code, name = "n") %>%
        left_join(dat0 %>% group_by(ancestry_code) %>%
                    summarise(families = n_distinct(FID), .groups = "drop"),
                  by = "ancestry_code"))

# -----------------------------
# Fit, with cluster-robust SEs
# -----------------------------
fit_one_ancestry <- function(df_anc, anc_label) {

  # standardise PRS WITHIN the analytic sample, after all filtering
  df_anc <- df_anc %>% mutate(across(all_of(prs_names), ~ as.numeric(scale(.x))))

  covar_rhs <- paste(c("sex", "batch", pc_names), collapse = " + ")

  fit_logistic <- function(outcome, prs) {
    fml <- as.formula(paste0(outcome, " ~ ", prs, " + ", covar_rhs))

    df2 <- df_anc %>%
      filter(!is.na(.data[[outcome]]),
             !if_any(all_of(c(pc_names, prs)), ~ is.na(.x)),
             !is.na(FID))

    if (nrow(df2) < 50) return(NULL)
    n_cases    <- sum(df2[[outcome]] == 1, na.rm = TRUE)
    n_controls <- sum(df2[[outcome]] == 0, na.rm = TRUE)
    if (n_cases == 0 || n_controls == 0) return(NULL)

    n_eff <- 4 / ((1 / n_cases) + (1 / n_controls))
    n_clusters <- n_distinct(df2$FID)

    fit <- tryCatch(glm(fml, data = df2, family = binomial()),
                    error = function(e) NULL)
    if (is.null(fit) || !isTRUE(fit$converged)) return(NULL)

    sm <- summary(fit)$coefficients
    if (!(prs %in% rownames(sm))) return(NULL)
    beta      <- sm[prs, "Estimate"]
    se_model  <- sm[prs, "Std. Error"]
    p_model   <- sm[prs, "Pr(>|z|)"]

    if (CLUSTER_SE) {
      ct <- tryCatch(
        lmtest::coeftest(fit,
          vcov = sandwich::vcovCL(fit, cluster = df2$FID, type = CLUSTER_TYPE)),
        error = function(e) NULL)
      if (is.null(ct) || !(prs %in% rownames(ct))) {
        # Do NOT fall back to unclustered model SEs — the analysis specifies
        # family-clustered SEs, so a model whose clustering fails is excluded.
        message("    EXCLUDED (clustered vcov failed): ", outcome, " ~ ", prs)
        return(NULL)
      } else {
        se <- ct[prs, "Std. Error"]; p <- ct[prs, "Pr(>|z|)"]
        se_source <- paste0("cluster-robust (", CLUSTER_TYPE, ")")
      }
    } else {
      se <- se_model; p <- p_model; se_source <- "model"
    }

    tibble(
      ancestry = anc_label, outcome = outcome, PRS = prs, model = "logistic",
      estimate = beta, std_error = se, p_value = p,
      std_error_model = se_model, p_value_model = p_model,
      se_inflation = se / se_model, se_source = se_source,
      OR = exp(beta),
      CI_low  = exp(beta - 1.96 * se),
      CI_high = exp(beta + 1.96 * se),
      n = nobs(fit), n_clusters = n_clusters,
      n_cases = n_cases, n_controls = n_controls, n_eff = n_eff)
  }

  map_dfr(outcomes_bin, function(o) map_dfr(prs_names, function(p) fit_logistic(o, p)))
}

# -----------------------------
# Run per ancestry
# -----------------------------
for (code in c("EUR", "AFR", "AMR")) {
  anc_label <- anc_map_full[[code]]
  df_anc <- dat0 %>% filter(ancestry_code == code)
  if (nrow(df_anc) == 0) { message("Skipping ", anc_label, ": no rows"); next }

  message("\nRunning ", anc_label, ": n=", nrow(df_anc),
          " in ", n_distinct(df_anc$FID), " families")

  res_anc <- fit_one_ancestry(df_anc, anc_label)
  if (is.null(res_anc) || nrow(res_anc) == 0) { message("No results"); next }

  res_anc <- res_anc %>%
    mutate(trait_code = case_when(
             str_starts(PRS, "EUR_")  ~ str_split_fixed(PRS, "_", 3)[, 2],
             str_starts(PRS, "META_") ~ str_split_fixed(PRS, "_", 3)[, 2],
             TRUE                     ~ str_split_fixed(PRS, "_", 3)[, 1]),
           trait = recode(trait_code, !!!trait_map, .default = trait_code)) %>%
    select(-trait_code) %>%
    group_by(trait) %>%
    mutate(q_BH = p.adjust(p_value, "BH")) %>%
    ungroup()

  out_path <- file.path(OUT_DIR,
    paste0("PRS_outcome_", gsub(" ", "_", anc_label), ".tsv"))
  fwrite(res_anc, out_path, sep = "\t")
  message("  wrote ", out_path,
          " | median SE inflation from clustering: ",
          sprintf("%.3f", median(res_anc$se_inflation, na.rm = TRUE)))
}

# =========================================
# Paper-ready tables
# =========================================
p_fmt <- function(x) ifelse(is.na(x), "",
            ifelse(x < 1e-6, formatC(x, format="e", digits=2),
              ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))))
q_fmt <- p_fmt
num2  <- function(x) ifelse(is.na(x), "", sprintf("%.2f", x))

read_ancestry_tsv <- function(anc_label) {
  f <- file.path(OUT_DIR, paste0("PRS_outcome_", gsub(" ", "_", anc_label), ".tsv"))
  if (!file.exists(f)) { message("No file for ", anc_label); return(NULL) }
  read_tsv(f, show_col_types = FALSE)
}

write_paper_table <- function(anc_label, filename) {
  df <- read_ancestry_tsv(anc_label)
  if (is.null(df) || !nrow(df)) return(invisible(NULL))

  out <- df %>%
    filter(
      !is.na(estimate), is.finite(estimate), abs(estimate) <= 10,
      !is.na(std_error), is.finite(std_error), std_error > 0, std_error <= 10
    ) %>%
    mutate(outcome_label = recode(outcome, !!!outcome_map, .default = outcome)) %>%
    arrange(trait, outcome_label) %>%
    mutate(
      `Effect (OR)` = paste0(num2(OR), " [", num2(CI_low), "-", num2(CI_high), "]"),
      Display = paste0("b ", num2(estimate), " (", num2(std_error), "); ",
                       "OR ", num2(OR), " [", num2(CI_low), "-", num2(CI_high), "]; ",
                       "p = ", p_fmt(p_value), "; BH q = ", q_fmt(q_BH),
                       ifelse(!is.na(q_BH) & q_BH < 0.05, " *", "")),
      `Significant (BH q<0.05)` = !is.na(q_BH) & q_BH < 0.05) %>%
    select(Ancestry = ancestry, Trait = trait, Outcome = outcome_label,
           `Effect (OR)`, Beta = estimate,
           `Standard error (clustered)` = std_error,
           `Standard error (model-based)` = std_error_model,
           `SE inflation` = se_inflation,
           `CI lower` = CI_low, `CI upper` = CI_high,
           `p-value` = p_value, `BH q-value` = q_BH, `Significant (BH q<0.05)`,
           `Sample size (n)` = n, `Families (n)` = n_clusters,
           Cases = n_cases, Controls = n_controls, Display)

  write_csv(out, filename)
  message("Wrote: ", filename)
  out
}

for (a in c("European", "African", "Admixed American")) {
  write_paper_table(a, file.path(OUT_DIR,
    paste0("prs_perinatal_clustered_", gsub(" ", "_", a), ".csv")))
}

# =====================================================
# Stouffer meta-analysis (N_eff weighted)
# =====================================================
files <- c(European = file.path(OUT_DIR, "PRS_outcome_European.tsv"),
           African  = file.path(OUT_DIR, "PRS_outcome_African.tsv"),
           `Admixed American` = file.path(OUT_DIR, "PRS_outcome_Admixed_American.tsv"))
files <- files[file.exists(files)]

res_all <- imap_dfr(files, function(path, anc) {
  read_tsv(path, show_col_types = FALSE) %>% mutate(ancestry = anc)
})

res_all_clean <- res_all %>%
  filter(trait %in% traits_keep) %>%
  mutate(outcome_label = recode(outcome, !!!outcome_map, .default = outcome),
         z_ancestry = estimate / std_error,      # uses the CLUSTERED SE
         weight_stouffer = sqrt(n_eff)) %>%
  filter(!is.na(estimate), is.finite(estimate), abs(estimate) <= 10,
         !is.na(std_error), is.finite(std_error), std_error > 0, std_error <= 10,
         is.finite(z_ancestry), is.finite(weight_stouffer), weight_stouffer > 0,
         n_cases > 0, n_controls > 0) %>%
  group_by(trait, outcome_label) %>%
  filter(n_distinct(ancestry) == 3) %>%
  ungroup()

write_csv(res_all_clean, file.path(OUT_DIR, "stouffer_meta_input.csv"))

meta_stouffer <- res_all_clean %>%
  group_by(trait, outcome_label) %>%
  summarise(
    z_meta = sum(weight_stouffer * z_ancestry, na.rm = TRUE) /
             sqrt(sum(weight_stouffer^2, na.rm = TRUE)),
    p_value_meta = 2 * pnorm(abs(z_meta), lower.tail = FALSE),
    n_total = sum(n, na.rm = TRUE),
    n_families_total = sum(n_clusters, na.rm = TRUE),
    n_cases_total = sum(n_cases, na.rm = TRUE),
    n_controls_total = sum(n_controls, na.rm = TRUE),
    n_eff_total = sum(n_eff, na.rm = TRUE),
    k_ancestries = n_distinct(ancestry),
    ancestries_included = paste(sort(unique(ancestry)), collapse = "; "),
    beta_neff_weighted = sum(n_eff * estimate, na.rm = TRUE) / sum(n_eff, na.rm = TRUE),
    OR_neff_weighted = exp(beta_neff_weighted),
    .groups = "drop") %>%
  group_by(trait) %>%
  mutate(q_BH_meta = p.adjust(p_value_meta, method = "BH")) %>%
  ungroup() %>%
  arrange(trait, outcome_label)

write_csv(meta_stouffer, file.path(OUT_DIR, "stouffer_meta.csv"))

meta_paper <- meta_stouffer %>%
  arrange(trait, outcome_label) %>%
  mutate(`Effect (OR)` = paste0(num2(OR_neff_weighted), " (descriptive; no CI)"),
         Display = paste0("b = ", num2(beta_neff_weighted), "; OR ",
                          num2(OR_neff_weighted), " (descriptive); Z = ", num2(z_meta),
                          "; p = ", p_fmt(p_value_meta), "; BH q = ", q_fmt(q_BH_meta),
                          ifelse(!is.na(q_BH_meta) & q_BH_meta < 0.05, " *", "")),
         `Significant (BH q<0.05)` = !is.na(q_BH_meta) & q_BH_meta < 0.05) %>%
  select(Trait = trait, Outcome = outcome_label,
         `Meta Z` = z_meta, `Meta p-value` = p_value_meta,
         `Meta BH q-value` = q_BH_meta, `Significant (BH q<0.05)`,
         `Beta, N_eff-weighted` = beta_neff_weighted,
         `OR, N_eff-weighted` = OR_neff_weighted, `Effect (OR)`,
         `N total` = n_total, `Families total` = n_families_total,
         `Cases total` = n_cases_total, `Controls total` = n_controls_total,
         `Effective sample size total` = n_eff_total,
         `Number of ancestries` = k_ancestries,
         `Ancestries included` = ancestries_included, Display)

write_csv(meta_paper, file.path(OUT_DIR, "stouffer_meta_paper_ready.csv"))

# -----------------------------
# Summary
# -----------------------------
cat("\n===== Sample: unaffected siblings with PRS | CLUSTER_SE:", CLUSTER_SE, "=====\n")
cat("\nEffect of clustering on standard errors (ratio clustered/model-based):\n")
print(as.data.frame(
  res_all %>% group_by(ancestry) %>%
    summarise(models = n(),
              median_inflation = median(se_inflation, na.rm = TRUE),
              max_inflation = max(se_inflation, na.rm = TRUE),
              sig_model = sum(p.adjust(p_value_model, "BH") < 0.05, na.rm = TRUE),
              sig_clustered = sum(q_BH < 0.05, na.rm = TRUE),
              .groups = "drop")))

cat("\nMeta-analysis, FDR-significant:\n")
print(as.data.frame(meta_stouffer %>% filter(q_BH_meta < 0.05) %>%
        select(trait, outcome_label, OR_neff_weighted, z_meta, q_BH_meta)))

message("\nWrote results to: ", OUT_DIR)
