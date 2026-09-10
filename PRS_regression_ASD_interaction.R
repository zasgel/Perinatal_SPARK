# =====================================================
# PRS ~ Perinatal Outcomes — POPULATION-LEVEL, PRS x ASD INTERACTION — SPARK V20
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

# -----------------------------
# Toggles
# -----------------------------
CLUSTER_TYPE <- "HC1"  
MIN_EVENTS   <- 5       # minimum observations required in each outcome x ASD cell

OUT_DIR <- file.path(PC_DIR, "population_regression_ASD_interaction_V20")
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

term_order <- c(
  "PRS (unaffected, ASD=0)",
  "ASD status",
  "PRS x ASD interaction"
)

# -----------------------------
# Phenotypes
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
    is.na(x) | x == ""     ~ 0L,          # not endorsed = control
    x == "1"               ~ 1L,          # endorsed = case
    x == "na_survey_logic" ~ NA_integer_, # not administered
    TRUE                   ~ NA_integer_
  )
}

basic_med <- fread(
  BASIC_MED_PATH,
  sep = ",",
  header = TRUE,
  quote = "\"",
  colClasses = "character"
) %>%
  as_tibble() %>%
  select(subject_sp_id, all_of(outcomes_bin)) %>%
  rename(IID = subject_sp_id) %>%
  mutate(
    IID = as.character(IID),
    across(all_of(outcomes_bin), clean_v20_binary)
  ) %>%
  distinct(IID, .keep_all = TRUE)

# -----------------------------
# V20 covariates + explicit child sample definition
# -----------------------------
blank_na <- function(x) {
  x <- trimws(as.character(x))
  dplyr::if_else(x %in% c("", "NA", "0", "."), NA_character_, x)
}

parse_asd <- function(x) {
  x <- tolower(trimws(as.character(x)))
  case_when(
    x %in% c("true","1","yes","y")  ~ 1L,
    x %in% c("false","0","no","n")  ~ 0L,
    TRUE ~ NA_integer_
  )
}

split_ids <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", ".")] <- NA_character_

  vals <- x[!is.na(x)]
  if (!length(vals)) return(character(0))

  ids <- unlist(strsplit(vals, "\\|"))
  ids <- trimws(ids)
  unique(ids[nzchar(ids)])
}

roles <- fread(ROLES_PATH, sep = ",", header = TRUE, quote = "\"") %>%
  as_tibble() %>%
  transmute(
    IID = as.character(subject_sp_id),
    biomother = blank_na(biomother_sp_id),
    biofather = blank_na(biofather_sp_id),
    affected_sibs = blank_na(affected_sibling_sp_id),
    control_sibs  = blank_na(control_sibling_sp_id)
  )

role_subject_ids      <- unique(roles$IID)
affected_sibling_ids  <- split_ids(roles$affected_sibs)
control_sibling_ids   <- split_ids(roles$control_sibs)

child_ids <- unique(c(
  role_subject_ids,
  affected_sibling_ids,
  control_sibling_ids
))

parent_ids <- unique(c(roles$biomother, roles$biofather))
parent_ids <- parent_ids[!is.na(parent_ids)]

indiv <- fread(INDIV_PATH, sep = ",", header = TRUE, quote = "\"") %>%
  as_tibble() %>%
  transmute(
    IID = as.character(subject_sp_id),
    FID = as.character(family_sf_id),
    sex = factor(sex),
    asd = parse_asd(asd),
    age_at_registration_months =
      suppressWarnings(as.numeric(age_at_registration_months))
  ) %>%
  distinct(IID, .keep_all = TRUE)


sibling_ids <- unique(c(affected_sibling_ids, control_sibling_ids))
missing_sibling_ids <- setdiff(sibling_ids, indiv$IID)
parent_child_overlap <- intersect(child_ids, parent_ids)

message("V20 explicit child-ID construction:")
message("  roles-row IDs: ", length(role_subject_ids))
message("  affected sibling IDs: ", length(affected_sibling_ids))
message("  control sibling IDs: ", length(control_sibling_ids))
message("  unique child IDs total: ", length(child_ids))
message("  sibling IDs absent from individuals_registration: ",
        length(missing_sibling_ids))
message("  explicit child IDs also listed as biological parents: ",
        length(parent_child_overlap))

if (length(missing_sibling_ids) > 0) {
  warning("Some sibling IDs listed in roles are absent from individuals_registration. Example(s): ",
          paste(head(missing_sibling_ids, 10), collapse = ", "))
}

meta_cov <- indiv %>%
  filter(IID %in% child_ids, !is.na(FID), !is.na(asd)) %>%
  distinct(IID, .keep_all = TRUE)

message("V20 analytic child sample with known ASD status: ", nrow(meta_cov),
        " individuals in ", n_distinct(meta_cov$FID), " families | ASD cases: ",
        sum(meta_cov$asd == 1, na.rm = TRUE),
        " | unaffected: ", sum(meta_cov$asd == 0, na.rm = TRUE))

print(meta_cov %>% count(asd, name = "n"))

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

IWES_META_PATH <- "/Users/asgelz01/Downloads/SPARK.iWES_v3.2024_08.sample_metadata.tsv"
batch_cov <- fread(IWES_META_PATH, sep = "\t", header = TRUE) %>%
  as_tibble() %>%
  transmute(IID = as.character(spid), batch = factor(batch)) %>%
  distinct(IID, .keep_all = TRUE)

# -----------------------------
# Merge
# -----------------------------
dat0 <- basic_med %>%
  inner_join(meta_cov, by = "IID") %>%
  inner_join(batch_cov, by = "IID") %>%
  inner_join(pcs, by = "IID") %>%
  inner_join(prs_wide, by = "IID")

prs_names <- setdiff(names(dat0),
  c("IID","FID", outcomes_bin, "sex","batch","asd","age_at_registration_months", "ancestry_code", pc_names))

message("\nAnalytic sample: ", nrow(dat0), " individuals in ",
        n_distinct(dat0$FID), " families | ASD cases: ",
        sum(dat0$asd == 1, na.rm = TRUE), " | unaffected: ",
        sum(dat0$asd == 0, na.rm = TRUE), " | unknown ASD status: ",
        sum(is.na(dat0$asd)))
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
    # Main effects + PRS x ASD interaction
    fml <- as.formula(
      paste0(outcome, " ~ ", prs, " * asd + sex + batch + ",
             paste(pc_names, collapse = " + "))
    )


    df2 <- df_anc %>%
      filter(!is.na(.data[[outcome]]),
             !if_any(all_of(c(pc_names, prs)), ~ is.na(.x)),
             !is.na(asd), !is.na(sex), !is.na(batch), !is.na(FID))

    if (nrow(df2) < 50) return(NULL)
    n_cases    <- sum(df2[[outcome]] == 1, na.rm = TRUE)
    n_controls <- sum(df2[[outcome]] == 0, na.rm = TRUE)


    cells <- table(
      factor(df2[[outcome]], levels = c(0, 1)),
      factor(df2$asd,        levels = c(0, 1))
    )

    if (any(cells < MIN_EVENTS)) {
      message(
        "    EXCLUDED (sparse outcome x ASD cell): ", outcome, " ~ ", prs,
        " * ASD | cells: ",
        paste(
          sprintf(
            "y%s/asd%s=%d",
            rep(c(0, 1), times = 2),
            rep(c(0, 1), each = 2),
            as.vector(cells)
          ),
          collapse = ", "
        )
      )
      return(NULL)
    }

    # Interaction requires both ASD groups; redundant after the 2x2 cell rule,
    # but retained as a defensive check.
    if (n_distinct(df2$asd) < 2) return(NULL)

    n_eff <- 4 / ((1 / n_cases) + (1 / n_controls))
    n_clusters <- n_distinct(df2$FID)

    fit <- tryCatch(glm(fml, data = df2, family = binomial()),
                    error = function(e) NULL)
    if (is.null(fit) || !fit$converged) return(NULL)

    sm <- summary(fit)$coefficients

    int1 <- paste0(prs, ":asd")
    int2 <- paste0("asd:", prs)
    int_term <- if (int1 %in% rownames(sm)) int1 else if (int2 %in% rownames(sm)) int2 else NA_character_


    terms_wanted <- c(`PRS (unaffected, ASD=0)` = prs,
                      `ASD status` = "asd",
                      `PRS x ASD interaction` = int_term)
    terms_wanted <- terms_wanted[!is.na(terms_wanted)]


    ct <- tryCatch(
      lmtest::coeftest(
        fit,
        vcov = sandwich::vcovCL(fit, cluster = df2$FID, type = CLUSTER_TYPE)
      ),
      error = function(e) NULL
    )
    if (is.null(ct)) {
      message("    EXCLUDED (clustered vcov failed): ", outcome, " ~ ", prs, " * ASD")
      return(NULL)
    }

    out_terms <- purrr::imap_dfr(terms_wanted, function(term, label) {
      if (!(term %in% rownames(ct)) || !(term %in% rownames(sm))) return(NULL)

      beta     <- ct[term, "Estimate"]
      se       <- ct[term, "Std. Error"]
      p        <- ct[term, "Pr(>|z|)"]
      se_model <- sm[term, "Std. Error"]
      p_model  <- sm[term, "Pr(>|z|)"]

      tibble(
        ancestry = anc_label,
        outcome = outcome,
        PRS = prs,
        term = label,
        model_term = term,
        model = "logistic_PRSxASD",
        estimate = beta,
        std_error = se,
        p_value = p,
        std_error_model = se_model,
        p_value_model = p_model,
        se_inflation = se / se_model,
        se_source = paste0("cluster-robust (", CLUSTER_TYPE, ")"),
        OR = exp(beta),
        CI_low = exp(beta - 1.96 * se),
        CI_high = exp(beta + 1.96 * se),
        n = nobs(fit),
        n_clusters = n_clusters,
        n_cases = n_cases,
        n_controls = n_controls,
        n_eff = n_eff,
        n_asd = sum(df2$asd == 1, na.rm = TRUE),
        n_unaffected = sum(df2$asd == 0, na.rm = TRUE)
      )
    })

    out_terms
  }

  map_dfr(outcomes_bin, function(o)
    map_dfr(prs_names, function(p) fit_logistic(o, p)))
}

# -----------------------------
# Run per ancestry
# -----------------------------
all_ancestry_results <- list()

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
    group_by(trait, term) %>%
    mutate(q_BH = p.adjust(p_value, "BH")) %>%
    ungroup() %>%
    mutate(
      outcome_label_sort = recode(outcome, !!!outcome_map, .default = outcome),
      term = factor(term, levels = term_order)
    ) %>%
    arrange(trait, outcome_label_sort, term) %>%
    mutate(term = as.character(term)) %>%
    select(-outcome_label_sort)

  all_ancestry_results[[anc_label]] <- res_anc

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
    mutate(
      outcome_label = recode(outcome, !!!outcome_map, .default = outcome),
      term = factor(term, levels = term_order)
    ) %>%
    arrange(trait, outcome_label, term) %>%
    mutate(
      term = as.character(term),
      `Effect (OR)` = paste0(num2(OR), " [", num2(CI_low), "-", num2(CI_high), "]"),
      Display = paste0("b ", num2(estimate), " (", num2(std_error), "); ",
                       "OR ", num2(OR), " [", num2(CI_low), "-", num2(CI_high), "]; ",
                       "p = ", p_fmt(p_value), "; BH q = ", q_fmt(q_BH),
                       ifelse(!is.na(q_BH) & q_BH < 0.05, " *", "")),
      `Significant (BH q<0.05)` = !is.na(q_BH) & q_BH < 0.05) %>%
    select(Ancestry = ancestry, Trait = trait, Outcome = outcome_label, Term = term,
           `Effect (OR)`, Beta = estimate,
           `Standard error (clustered)` = std_error,
           `Standard error (model-based)` = std_error_model,
           `SE inflation` = se_inflation,
           `CI lower` = CI_low, `CI upper` = CI_high,
           `p-value` = p_value, `BH q-value` = q_BH, `Significant (BH q<0.05)`,
           `Sample size (n)` = n, `Families (n)` = n_clusters,
           ASD = n_asd, Unaffected = n_unaffected, Cases = n_cases, Controls = n_controls, Display)

  write_csv(out, filename)
  message("Wrote: ", filename)
  out
}

for (a in c("European", "African", "Admixed American")) {
  write_paper_table(a, file.path(OUT_DIR,
    paste0("prs_perinatal_clustered_", gsub(" ", "_", a), ".csv")))
}

# -----------------------------
# Summary
# -----------------------------
res_all <- bind_rows(all_ancestry_results)

if (nrow(res_all) == 0) {
  stop("No ancestry-specific models were successfully fitted; no summary can be produced.")
}

cat("\n===== Sample: probands + siblings | PRS x ASD interaction |",
    "cluster-robust SEs (", CLUSTER_TYPE, ") on family ID =====\n")
cat("\nEffect of clustering on standard errors (ratio clustered/model-based):\n")
print(as.data.frame(
  res_all %>% group_by(ancestry) %>%
    summarise(models = n(),
              median_inflation = median(se_inflation, na.rm = TRUE),
              max_inflation = max(se_inflation, na.rm = TRUE),
              sig_model = sum(p.adjust(p_value_model, "BH") < 0.05, na.rm = TRUE),
              sig_clustered = sum(q_BH < 0.05, na.rm = TRUE),
              .groups = "drop")))

cat("\nFDR-significant results (BH q < 0.05), by ancestry and term:\n")
sig <- res_all %>%
  filter(!is.na(q_BH), q_BH < 0.05) %>%
  mutate(outcome_label = recode(outcome, !!!outcome_map, .default = outcome)) %>%
  select(ancestry, term, trait, outcome_label, OR, CI_low, CI_high, p_value, q_BH) %>%
  arrange(ancestry, trait, outcome_label, factor(term, levels = term_order))
if (nrow(sig)) print(as.data.frame(sig)) else cat("  none\n")

cat("\nInteraction terms, ranked by p-value (top 15):\n")
int_tb <- res_all %>%
  filter(term == "PRS x ASD interaction") %>%
  mutate(outcome_label = recode(outcome, !!!outcome_map, .default = outcome)) %>%
  select(ancestry, trait, outcome_label, n_asd, n_unaffected,
         estimate, std_error, p_value, q_BH) %>%
  arrange(p_value)
if (nrow(int_tb)) print(as.data.frame(head(int_tb, 15))) else cat("  none fitted\n")

message("\nWrote results to: ", OUT_DIR)
