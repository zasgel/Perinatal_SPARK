# =========================================================
# Within-siblings PRS analysis — Per ancestry (EUR, AFR, AMR)
# + Stouffer meta-analysis across ancestries
# Separate meta for PRS_within and PRS_mean
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
BASIC_MED_PATH <- "/Users/asgelz01/Downloads/SPARKDataRelease_2025-07-14/basic_medical_screening-2025-07-14.csv"

OUT_DIR <- "/Users/asgelz01/Downloads/PRS_all/within_sibs_results_per_ancestry"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Toggles
# -----------------------------
KEEP_BETWEEN <- TRUE
MIN_EVENTS   <- 5
DIGITS_MZ    <- 6

# -----------------------------
# Shared labels
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

traits_keep <- c("Schizophrenia", "Obesity", "Asthma", "Major depression")

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

outcomes <- names(outcome_map)

anc_map_full <- c(EUR = "European", AFR = "African", AMR = "Admixed American")

# -----------------------------
# Load shared data
# -----------------------------
meta <- fread(META_PATH, sep = "\t") %>%
  rename(IID = spid, FID = sfid) %>%
  select(IID, FID, father, mother, sex, batch, identical_twins) %>%
  mutate(
    identical_twins = !(identical_twins %in% c(NA, "", "False", "0", "NA")),
    sex   = factor(sex),
    batch = factor(batch)
  )

basic_med <- fread(BASIC_MED_PATH) %>%
  select(subject_sp_id, all_of(outcomes)) %>%
  rename(IID = subject_sp_id) %>%
  mutate(
    across(all_of(outcomes), ~ tidyr::replace_na(., 0L)),
    across(all_of(outcomes), as.integer)
  )


# -----------------------------
# Formatting helpers for paper-ready tables 
# -----------------------------
num2  <- function(x, digits = 2, na = "") {
  ifelse(is.na(x), na, formatC(x, format = "f", digits = digits))
}
p_fmt <- function(p, na = "") {
  ifelse(is.na(p), na,
         ifelse(p < 1e-4,
                formatC(p, format = "e", digits = 2),
                formatC(p, format = "f", digits = 3)))
}
q_fmt <- p_fmt


build_per_ancestry_paper_table <- function(res_anc, term_keep, term_label) {
  res_anc %>%
    dplyr::filter(term == term_keep | is.na(term)) %>%
    dplyr::mutate(
      term = ifelse(is.na(term), term_keep, term),
      sig_bh = !is.na(q_BH) & q_BH < 0.05
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      `Effect (OR)` = ifelse(
        is.na(OR), "NA",
        paste0(num2(OR), " [", num2(CI_low), "–", num2(CI_high), "]")
      ),
      Display = paste0(
        "β ", ifelse(is.na(estimate), "NA", num2(estimate)), " (",
        ifelse(is.na(std.error), "NA", num2(std.error)), "); ",
        "OR ", ifelse(is.na(OR), "NA", num2(OR)),
        ifelse(is.na(OR), "", paste0(" [", num2(CI_low), "–", num2(CI_high), "]")),
        "; p = ", ifelse(is.na(p), "NA", p_fmt(p)),
        "; BH q = ", ifelse(is.na(q_BH), "NA", q_fmt(q_BH)),
        ifelse(!is.na(q_BH) & q_BH < 0.05, " *", "")
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Term = term_label) %>%
    dplyr::arrange(trait_full, outcome_label) %>%
    dplyr::select(
      `Trait`             = trait_full,
      `Outcome`           = outcome_label,
      `Term`,
      `Effect (OR)`,
      `Beta`              = estimate,
      `Standard error`    = std.error,
      `Odds ratio`        = OR,
      `CI lower`          = CI_low,
      `CI upper`          = CI_high,
      `p-value`           = p,
      `BH q-value`        = q_BH,
      `BH significant`    = sig_bh,
      `Individuals (n)`   = n_ind,
      `Families (n)`      = n_fam,
      `Events`            = events,
      `Non-events`        = nonevents,
      `Model note`        = note,
      Display
    )
}

# -----------------------------
# Helpers 
# -----------------------------
read_sscore_simple <- function(path) {
  df <- data.table::fread(path, sep = "\t", header = TRUE, check.names = FALSE)
  data.table::setnames(df, sub("^#", "", names(df)))
  if (!("IID" %in% names(df))) stop("IID missing in ", path)
  score_col <- if ("SCORESUM" %in% names(df)) "SCORESUM" else
    if ("SCORE1_SUM" %in% names(df)) "SCORE1_SUM" else
      names(df)[ncol(df)]
  tibble(IID = df$IID, PRS_raw = as.numeric(df[[score_col]]))
}

keep_one_per_mz_pair <- function(df, digits = DIGITS_MZ) {
  if (!nrow(df)) return(df)
  split_list <- split(df, df$FID, drop = TRUE)
  processed <- lapply(split_list, function(d) {
    if (!any(isTRUE(d$identical_twins))) return(d)
    if (!("PRS_z" %in% names(d))) return(d)
    tw <- d[isTRUE(d$identical_twins), , drop = FALSE]
    if (!nrow(tw)) return(d)
    grp <- round(tw$PRS_z, digits = digits)
    to_drop <- character(0)
    for (g in unique(grp)) {
      ids <- tw$IID[grp == g]
      if (length(ids) > 1) to_drop <- c(to_drop, ids[-1])
    }
    d[!(d$IID %in% to_drop), , drop = FALSE]
  })
  dplyr::bind_rows(processed)
}

prep_within_family_data <- function(prs_tbl, meta_anc, basic_med) {
  inner_join(meta_anc, basic_med, by = "IID") %>%
    inner_join(prs_tbl, by = "IID") %>%
    filter(!is.na(FID), !is.na(sex), !is.na(batch)) %>%
    add_count(FID, name = "fam_n") %>% filter(fam_n >= 2) %>% select(-fam_n) %>%
    mutate(PRS_z = as.numeric(scale(PRS_raw))) %>%
    group_by(FID) %>%
    mutate(
      PRS_mean   = mean(PRS_z, na.rm = TRUE),
      PRS_within = PRS_z - PRS_mean
    ) %>%
    ungroup() %>%
    keep_one_per_mz_pair() %>%
    add_count(FID, name = "fam_n") %>% filter(fam_n >= 2) %>% select(-fam_n)
}

get_safe_vcov <- function(fit) {
  V <- tryCatch(as.matrix(vcov(fit, use.hessian = FALSE)), error = function(e) NULL)
  if (!is.null(V) && all(is.finite(diag(V))) && all(diag(V) >= 0)) return(V)
  V <- tryCatch(as.matrix(vcov(fit)), error = function(e) NULL)
  if (!is.null(V) && all(is.finite(diag(V))) && all(diag(V) >= 0)) return(V)
  NULL
}

fit_one_condition <- function(yvar, d) {
  empty <- function(note) tibble(
    condition = yvar, term = NA_character_, estimate = NA_real_, std.error = NA_real_,
    z = NA_real_, p = NA_real_, OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
    n_ind = if (is.null(d)) 0L else nrow(d),
    n_fam = if (is.null(d)) 0L else dplyr::n_distinct(d$FID),
    events = NA_integer_, nonevents = NA_integer_, note = note
  )
  
  if (is.null(d)) return(empty("No data"))
  
  dd <- d %>% filter(!is.na(.data[[yvar]]))
  if (nrow(dd) < 4) return(empty("Too few rows"))
  
  dd <- dd %>% add_count(FID, name = "fam_n") %>% filter(fam_n >= 2) %>% select(-fam_n)
  if (nrow(dd) < 4) return(empty("No families with >=2 sibs"))
  
  fam_var <- dd %>% group_by(FID) %>%
    summarise(var_y = var(.data[[yvar]], na.rm = TRUE), .groups = "drop")
  if (all(is.na(fam_var$var_y) | fam_var$var_y == 0))
    return(empty("No within-family variation in outcome"))
  
  events    <- sum(dd[[yvar]] == 1, na.rm = TRUE)
  nonevents <- sum(dd[[yvar]] == 0, na.rm = TRUE)
  if (!is.na(events) && events <= MIN_EVENTS)
    return(empty(paste0("Too few events (<= ", MIN_EVENTS, ")")))
  
  # Separation guard
  form <- as.formula(paste0(yvar, " ~ PRS_mean + PRS_within + sex + batch + (1|FID)"))
  fit <- tryCatch(
    glmer(form, data = dd, family = binomial(),
          control = glmerControl(
            optimizer = "bobyqa", optCtrl = list(maxfun = 2e5),
            calc.derivs = TRUE,
            check.conv.grad = "ignore",
            check.conv.singular = "ignore",
            check.conv.hess = "ignore"
          ), nAGQ = 0),
    error = function(e) e
  )
  
  n_ind <- nrow(dd); n_fam <- dplyr::n_distinct(dd$FID)
  if (inherits(fit, "error"))
    return(tibble(condition=yvar, term=NA_character_, estimate=NA_real_, std.error=NA_real_,
                  z=NA_real_, p=NA_real_, OR=NA_real_, CI_low=NA_real_, CI_high=NA_real_,
                  n_ind=n_ind, n_fam=n_fam, events=events, nonevents=nonevents,
                  note=paste("Model error:", conditionMessage(fit))))
  
  fe <- lme4::fixef(fit)
  
  # Separation guard: drop implausible estimates
  if (any(abs(fe) > 10, na.rm = TRUE))
    return(tibble(condition=yvar, term=NA_character_, estimate=NA_real_, std.error=NA_real_,
                  z=NA_real_, p=NA_real_, OR=NA_real_, CI_low=NA_real_, CI_high=NA_real_,
                  n_ind=n_ind, n_fam=n_fam, events=events, nonevents=nonevents,
                  note="Possible separation: |beta| > 10"))
  
  terms <- names(fe)
  V <- get_safe_vcov(fit)
  
  if (!is.null(V)) {
    se      <- sqrt(pmax(0, diag(V)))
    z_val   <- ifelse(se > 0, fe / se, NA_real_)
    p_val   <- ifelse(is.finite(z_val), 2 * pnorm(abs(z_val), lower.tail = FALSE), NA_real_)
    OR      <- exp(fe); CI_low <- exp(fe - 1.96*se); CI_high <- exp(fe + 1.96*se)
    note_msg <- NA_character_
  } else {
    se <- z_val <- p_val <- OR <- CI_low <- CI_high <- rep(NA_real_, length(fe))
    note_msg <- "Non-PD vcov"
  }
  
  tibble(condition=yvar, term=terms, estimate=as.numeric(fe),
         std.error=as.numeric(se), z=as.numeric(z_val), p=as.numeric(p_val),
         OR=as.numeric(OR), CI_low=as.numeric(CI_low), CI_high=as.numeric(CI_high),
         n_ind=n_ind, n_fam=n_fam, events=events, nonevents=nonevents, note=note_msg) %>%
    filter(term %in% c("PRS_within", "PRS_mean"))
}


files_all <- list.files(SSCORE_DIR, pattern = "\\.sscore$", full.names = TRUE)
files_all <- files_all[!grepl("^META", basename(files_all), ignore.case = TRUE)]

label_file <- function(fp) {
  b <- tools::file_path_sans_ext(basename(fp))
  parts <- str_split_fixed(b, "_", 2)  # only 2 parts now, no _ASD suffix
  
  p1 <- parts[1]
  p2 <- parts[2]
  
  anc_codes   <- c("EUR", "AFR", "AMR")
  trait_codes <- names(trait_map)  # AS, OB, SCZ, MDD, PPH, PE, GD
  
  if (p1 %in% anc_codes && p2 %in% trait_codes) {
    # EUR_SCZ pattern
    return(tibble(ancestry = p1, trait = p2, file = fp))
  }
  if (p1 %in% trait_codes && p2 %in% anc_codes) {
    # SCZ_AFR pattern
    return(tibble(ancestry = p2, trait = p1, file = fp))
  }
  return(NULL)
}

labels <- map_dfr(files_all, label_file)

if (nrow(labels) == 0) stop("No PRS files matched expected naming patterns.")

# -----------------------------
# Run per ancestry
# -----------------------------
all_anc_res <- list()

for (code in c("EUR", "AFR", "AMR")) {
  anc_label  <- anc_map_full[[code]]
  labels_anc <- labels %>% filter(ancestry == code)
  # No PC-based ancestry filter on meta: the .sscore file is already
  # ancestry-restricted upstream, so the inner-join with prs_tbl defines
  # the per-ancestry cohort.
  meta_anc   <- meta

  if (nrow(labels_anc) == 0) {
    message("No PRS files for ancestry: ", code); next
  }
  if (nrow(meta_anc) == 0) {
    message("No individuals in meta (check META_PATH): ", code); next
  }

  message("\n=== Ancestry: ", anc_label,
          " | n_PRS_files=", nrow(labels_anc), " ===")
  
  traits_anc <- unique(labels_anc$trait)
  res_anc_list <- list()
  
  for (tr in traits_anc) {
    files_tr <- labels_anc$file[labels_anc$trait == tr]
    message(sprintf("  trait=%s | n_files=%d", tr, length(files_tr)))
    
    prs_tbl <- map_dfr(files_tr, read_sscore_simple) %>%
      distinct(IID, .keep_all = TRUE)
    
    d_tr <- prep_within_family_data(prs_tbl, meta_anc, basic_med)
    
    res_tr <- map_dfr(outcomes, ~ fit_one_condition(.x, d_tr)) %>%
      mutate(trait = tr, ancestry = anc_label)
    
    if (!KEEP_BETWEEN) res_tr <- res_tr %>% filter(term == "PRS_within")
    
    res_anc_list[[tr]] <- res_tr
  }
  
  res_anc <- bind_rows(res_anc_list) %>%
    mutate(
      trait_full    = recode(trait, !!!trait_map, .default = trait),
      outcome_label = recode(condition, !!!outcome_map, .default = condition)
    ) %>%
    group_by(trait_full, term) %>%
    mutate(q_BH = p.adjust(p, method = "BH")) %>%
    ungroup()
  
  out_path <- file.path(OUT_DIR, paste0("within_sibs_", gsub(" ", "_", anc_label), ".tsv"))
  fwrite(res_anc, out_path, sep = "\t")
  message("Wrote: ", out_path)


  term_labels <- c(
    PRS_within = "Within-family (PRS_within)",
    PRS_mean   = "Between-family (PRS_mean)"
  )
  for (term_keep in c("PRS_within", "PRS_mean")) {
    paper_anc <- build_per_ancestry_paper_table(
      res_anc, term_keep, term_labels[[term_keep]]
    )
    paper_tag <- if (term_keep == "PRS_within") "within" else "between"
    paper_path <- file.path(
      OUT_DIR,
      paste0("paper_ready_", paper_tag, "_sibs_",
             gsub(" ", "_", anc_label), ".csv")
    )
    readr::write_csv(paper_anc, paper_path)
    message("Wrote paper-ready: ", paper_path)
  }

  all_anc_res[[anc_label]] <- res_anc
}

# =====================================================
# Stouffer meta-analysis — within-family (PRS_within)
# and between-family (PRS_mean) separately
# =====================================================

# Combine all ancestry results
res_all_terms <- imap_dfr(all_anc_res, function(df, anc_name) {
  df %>% mutate(ancestry = anc_name)
}) %>%
  filter(trait_full %in% traits_keep)

# Template of all trait × outcome combinations to report
meta_template <- res_all_terms %>%
  distinct(trait_full, outcome_label) %>%
  tidyr::crossing(term = c("PRS_within", "PRS_mean"))

# Valid rows only for Stouffer calculation
res_meta_input <- res_all_terms %>%
  filter(
    term %in% c("PRS_within", "PRS_mean"),
    !is.na(estimate),
    !is.na(z),
    !is.na(std.error),
    std.error > 0,
    !is.na(events),
    !is.na(nonevents),
    events > 0,
    nonevents > 0
  ) %>%
  mutate(
    n_eff = 4 / ((1 / events) + (1 / nonevents)),
    z_score = z
  ) %>%
  filter(
    !is.na(n_eff),
    n_eff > 0,
    !is.na(z_score)
  )

res_meta_input <- res_meta_input %>%
  group_by(trait_full, outcome_label, term) %>%
  filter(n_distinct(ancestry) == 3) %>%
  ungroup()

# Save valid meta input for transparency
fwrite(
  res_meta_input,
  file.path(OUT_DIR, "stouffer_within_sibs_meta_input.tsv"),
  sep = "\t"
)

# Run Stouffer separately for each term
run_stouffer <- function(df_term, term_label) {
  
  df_term %>%
    mutate(
      weight_stouffer = if (term_label == "PRS_within") sqrt(n_fam) else sqrt(n_eff),
      weight_beta     = if (term_label == "PRS_within") n_fam else n_eff
    ) %>%
    filter(
      !is.na(weight_stouffer), weight_stouffer > 0,
      !is.na(weight_beta),     weight_beta > 0
    ) %>%
    group_by(trait_full, outcome_label) %>%
    summarise(
      term = term_label,
      z_meta = sum(weight_stouffer * z_score, na.rm = TRUE) /
        sqrt(sum(weight_stouffer^2, na.rm = TRUE)),
      p_value_meta     = 2 * pnorm(abs(z_meta), lower.tail = FALSE),
      n_ind_total      = sum(n_ind, na.rm = TRUE),
      n_fam_total      = sum(n_fam, na.rm = TRUE),
      n_cases_total    = sum(events, na.rm = TRUE),
      n_controls_total = sum(nonevents, na.rm = TRUE),
      n_eff_total      = sum(n_eff, na.rm = TRUE),
      k_ancestries     = n_distinct(ancestry),
      ancestries_included = paste(sort(unique(ancestry)), collapse = "; "),
      beta_weighted    = sum(weight_beta * estimate, na.rm = TRUE) /
        sum(weight_beta, na.rm = TRUE),
      OR_weighted      = exp(beta_weighted),
      .groups = "drop"
    ) %>%
    group_by(trait_full) %>%
    mutate(q_BH_meta = p.adjust(p_value_meta, method = "BH")) %>%
    ungroup() %>%
    arrange(trait_full, outcome_label)
}

meta_within <- run_stouffer(
  res_meta_input %>% filter(term == "PRS_within"),
  "PRS_within"
)

meta_between <- run_stouffer(
  res_meta_input %>% filter(term == "PRS_mean"),
  "PRS_mean"
)

write_csv(
  meta_within,
  file.path(OUT_DIR, "stouffer_meta_within_family.csv")
)

write_csv(
  meta_between,
  file.path(OUT_DIR, "stouffer_meta_between_family.csv")
)

message("Wrote Stouffer meta results for trait/outcome combinations with all 3 ancestries available.")
print(meta_within)
print(meta_between)

# =========================================
# Paper-ready table for Stouffer meta-analysis
# =========================================
build_meta_paper_table <- function(meta_df, term_label) {
  # Weighting scheme differs by term (set in run_stouffer above)
  weight_descr <- if (term_label == "PRS_within") "n_fam-weighted" else "n_eff-weighted"
  
  meta_df %>%
    mutate(
      `Effect (OR)` = paste0(
        num2(OR_weighted),
        " (descriptive, ", weight_descr, "; no CI\u2020)"
      ),
      Display = paste0(
        "\u03b2 = ", num2(beta_weighted), "; ",
        "OR ", num2(OR_weighted), " (descriptive, ", weight_descr, "); ",
        "Z = ", num2(z_meta), "; ",
        "p = ", p_fmt(p_value_meta), "; ",
        "BH q = ", q_fmt(q_BH_meta),
        ifelse(!is.na(q_BH_meta) & q_BH_meta < 0.05, " *", "")
      ),
      `Significant (BH q<0.05)` = !is.na(q_BH_meta) & q_BH_meta < 0.05,
      Term = term_label
    ) %>%
    select(
      `Trait`                          = trait_full,
      `Outcome`                        = outcome_label,
      `Term`,
      
      `Meta Z`                         = z_meta,
      `Meta p-value`                   = p_value_meta,
      `Meta BH q-value`                = q_BH_meta,
      `Significant (BH q<0.05)`,
      
      `Beta, weighted`                 = beta_weighted,
      `OR, weighted`                   = OR_weighted,
      `Effect (OR)`,
      
      `N individuals total`            = n_ind_total,
      `N families total`               = n_fam_total,
      `Cases total`                    = n_cases_total,
      `Controls total`                 = n_controls_total,
      `Effective sample size total`    = n_eff_total,
      
      `Number of ancestries`           = k_ancestries,
      `Ancestries included`            = ancestries_included,
      
      Display
    )
}

# ---- Paper-ready Stouffer meta tables ----
meta_within_paper  <- build_meta_paper_table(meta_within,  "PRS_within")
meta_between_paper <- build_meta_paper_table(meta_between, "PRS_mean")

write_csv(
  meta_within_paper,
  file.path(OUT_DIR, "paper_ready_stouffer_meta_within_family.csv")
)
write_csv(
  meta_between_paper,
  file.path(OUT_DIR, "paper_ready_stouffer_meta_between_family.csv")
)

message("Wrote paper-ready Stouffer meta tables.")
print(meta_within_paper)
print(meta_between_paper)
