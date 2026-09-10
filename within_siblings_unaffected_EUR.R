# =========================================================
# Within-siblings PRS analysis — UNAFFECTED SIBLINGS ONLY — EUR
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

OUT_DIR <- "/Users/asgelz01/Downloads/PRS_all/within_sibs_results_unaffected_only_EUR"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Toggles
# -----------------------------
KEEP_BETWEEN <- TRUE
MIN_EVENTS   <- 5    
DIGITS_MZ    <- 6    

BH_SCOPE <- "within_trait"

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

anc_map_full <- c(EUR = "European")

# -----------------------------
# Load shared data
# -----------------------------
blank_na <- function(x) {
  x <- trimws(as.character(x))
  dplyr::if_else(x %in% c("", "NA", "0", "."), NA_character_, x)
}

meta_raw <- fread(META_PATH, sep = "\t") %>%
  as_tibble() %>%
  mutate(across(c(spid, sfid), as.character),
         father = blank_na(father), mother = blank_na(mother))

# -----------------------------
# MZ twin identification (from the identical_twins column itself)
# -----------------------------
parse_twin_tokens <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "0", ".", "False", "FALSE", "false",
             "True", "TRUE", "true", "1")] <- NA_character_
  strsplit(x, "[,;/|[:space:]]+")
}

build_mz_groups <- function(ids, twin_col) {
  toks <- parse_twin_tokens(twin_col)
  edges <- map2_dfr(ids, toks, function(i, t) {
    t <- t[!is.na(t) & nzchar(t) & t != i]
    if (!length(t)) return(NULL)
    tibble(a = i, b = t)
  })
  if (!nrow(edges)) return(tibble(IID = character(0), mz_group = character(0)))

  nodes  <- unique(c(edges$a, edges$b))
  parent <- setNames(nodes, nodes)
  find_root <- function(x) {
    while (parent[[x]] != x) {
      parent[[x]] <<- parent[[parent[[x]]]]  # path compression
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

n_nonblank_twin <- sum(!is.na(blank_na(meta_raw$identical_twins)) &
                       !(trimws(meta_raw$identical_twins) %in%
                           c("True","TRUE","true","False","FALSE","false","1")))
frac_known <- if (nrow(mz_groups)) mean(mz_groups$IID %in% meta_raw$spid) else 0

if (nrow(mz_groups) > 0 && frac_known > 0.5) {
  MZ_MODE <- "spid"
  message("MZ mode: SPID-linked. ", nrow(mz_groups), " individuals in ",
          n_distinct(mz_groups$mz_group), " MZ sets (",
          round(100 * frac_known, 1), "% of linked IDs found in metadata).")
} else {
  MZ_MODE <- "prs_fallback"
  message("MZ mode: FALLBACK (PRS similarity). identical_twins does not appear ",
          "to contain co-twin SPIDs (", n_nonblank_twin,
          " non-boolean non-blank values). Inspect the column before trusting this:")
  print(head(table(meta_raw$identical_twins, useNA = "ifany"), 10))
}

# -----------------------------
# V20 explicit proband + sibling construction
# -----------------------------
parse_asd <- function(x) {
  x <- tolower(trimws(as.character(x)))
  case_when(x %in% c("true","1","yes","y") ~ 1L,
            x %in% c("false","0","no","n") ~ 0L,
            TRUE ~ NA_integer_)
}

split_ids <- function(x) {
  x <- trimws(as.character(x)); x[x %in% c("", "NA", ".")] <- NA_character_
  vals <- x[!is.na(x)]; if (!length(vals)) return(character(0))
  ids <- trimws(unlist(strsplit(vals, "\\|")))
  unique(ids[nzchar(ids)])
}

roles <- fread(ROLES_PATH) %>% as_tibble() %>% transmute(
  IID = as.character(subject_sp_id),
  affected_sibs = blank_na(affected_sibling_sp_id),
  control_sibs = blank_na(control_sibling_sp_id))
affected_sibling_ids <- split_ids(roles$affected_sibs)
control_sibling_ids <- split_ids(roles$control_sibs)
child_ids <- unique(c(roles$IID, affected_sibling_ids, control_sibling_ids))

indiv <- fread(INDIV_PATH) %>% as_tibble() %>% transmute(
  IID = as.character(subject_sp_id), FID = as.character(family_sf_id),
  sex = factor(sex), asd = parse_asd(asd)) %>% distinct(IID, .keep_all = TRUE)
missing_sibling_ids <- setdiff(unique(c(affected_sibling_ids, control_sibling_ids)), indiv$IID)
message("Sibling IDs absent from individuals_registration: ", length(missing_sibling_ids))

genomic_cov <- meta_raw %>% transmute(
  IID = as.character(spid), batch = factor(batch),
  mz_flag = !(identical_twins %in% c(NA, "", "False", "FALSE", "false", "0", "NA"))) %>%
  left_join(mz_groups, by="IID") %>% distinct(IID, .keep_all=TRUE)

meta <- indiv %>% filter(IID %in% child_ids, !is.na(FID), !is.na(sex), !is.na(asd), asd == 0L) %>%
  inner_join(genomic_cov, by="IID") %>%
  select(IID, FID, sex, batch, asd, mz_flag, mz_group) %>% distinct(IID, .keep_all=TRUE)
message("V20 unaffected proband/sibling sample: ", nrow(meta), " individuals in ", n_distinct(meta$FID), " families.")

# V20 phenotype coding: blank=0; 1=case; na_survey_logic=missing/not administered
clean_v20_binary <- function(x) {
  x <- trimws(as.character(x))
  case_when(is.na(x) | x == "" ~ 0L, x == "1" ~ 1L, x == "na_survey_logic" ~ NA_integer_, TRUE ~ NA_integer_)
}
basic_med <- fread(BASIC_MED_PATH, colClasses="character") %>% as_tibble() %>%
  select(subject_sp_id, all_of(outcomes)) %>% rename(IID=subject_sp_id) %>%
  mutate(IID=as.character(IID), across(all_of(outcomes), clean_v20_binary)) %>%
  distinct(IID, .keep_all=TRUE)

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
        paste0(num2(OR), " [", num2(CI_low), "\u2013", num2(CI_high), "]")
      ),
      Display = paste0(
        "\u03b2 ", ifelse(is.na(estimate), "NA", num2(estimate)), " (",
        ifelse(is.na(std.error), "NA", num2(std.error)), "); ",
        "OR ", ifelse(is.na(OR), "NA", num2(OR)),
        ifelse(is.na(OR), "", paste0(" [", num2(CI_low), "\u2013", num2(CI_high), "]")),
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
  tibble(IID = as.character(df$IID), PRS_raw = as.numeric(df[[score_col]]))
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
    mutate(
      PRS_mean   = mean(PRS_z, na.rm = TRUE),
      PRS_within = PRS_z - PRS_mean
    ) %>%
    ungroup()
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

  # families where the outcome varies -- these identify PRS_within
  .fam <- dd %>% group_by(FID) %>%
    summarise(n1 = sum(.data[[yvar]] == 1L), n0 = sum(.data[[yvar]] == 0L),
              .groups = "drop") %>%
    mutate(disc = n1 > 0 & n0 > 0)
  n_disc_fam  <- sum(.fam$disc)
  events_disc <- sum(.fam$n1[.fam$disc])

  events    <- sum(dd[[yvar]] == 1, na.rm = TRUE)
  nonevents <- sum(dd[[yvar]] == 0, na.rm = TRUE)

  if (events < MIN_EVENTS || nonevents < MIN_EVENTS)
    return(empty(paste0("Too few events or non-events (< ", MIN_EVENTS, ")")))

  if (events_disc < MIN_EVENTS)
    return(empty(paste0("Too few events in discordant families (< ", MIN_EVENTS, ")")))

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


  prs_terms <- intersect(c("PRS_mean", "PRS_within"), names(fe))
  if (any(abs(fe[prs_terms]) > 10, na.rm = TRUE))
    return(tibble(condition=yvar, term=NA_character_, estimate=NA_real_, std.error=NA_real_,
                  z=NA_real_, p=NA_real_, OR=NA_real_, CI_low=NA_real_, CI_high=NA_real_,
                  n_ind=n_ind, n_fam=n_fam, events=events, nonevents=nonevents,
                  note="Possible separation/instability in PRS term: |beta| > 10"))

  terms <- names(fe)
  V <- get_safe_vcov(fit)

  if (!is.null(V)) {
    se <- sqrt(pmax(0, diag(V)))
    names(se) <- names(fe)
    prs_se <- se[prs_terms]
    if (any(!is.finite(prs_se)) || any(prs_se <= 0) || any(prs_se > 10)) {
      return(tibble(condition=yvar, term=NA_character_, estimate=NA_real_, std.error=NA_real_,
                    z=NA_real_, p=NA_real_, OR=NA_real_, CI_low=NA_real_, CI_high=NA_real_,
                    n_ind=n_ind, n_fam=n_fam, events=events, nonevents=nonevents,
                    note="Unstable PRS term: non-finite, non-positive, or SE > 10"))
    }
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
         n_ind=n_ind, n_fam=n_fam, n_disc_fam=n_disc_fam, events_disc=events_disc,
         events=events, nonevents=nonevents, note=note_msg) %>%
    filter(term %in% c("PRS_within", "PRS_mean"))
}

# -----------------------------
# BH correction (explicit, consistent testing family)
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


files_all <- list.files(SSCORE_DIR, pattern = "\\.sscore$", full.names = TRUE)

label_file <- function(fp) {
  b <- tools::file_path_sans_ext(basename(fp))
  parts <- str_split_fixed(b, "_", 2)  # only 2 parts now, no _ASD suffix

  p1 <- parts[1]
  p2 <- parts[2]

  anc_codes   <- c("EUR")
  trait_codes <- names(trait_map)  # AS, OB, SCZ, MDD, PPH, PE, GD

  if (p1 %in% anc_codes && p2 %in% trait_codes) {
    return(tibble(ancestry = p1, trait = p2, file = fp))
  }
  if (p1 %in% trait_codes && p2 %in% anc_codes) {
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

for (code in c("EUR")) {
  anc_label  <- anc_map_full[[code]]
  labels_anc <- labels %>% filter(ancestry == code)
  meta_anc   <- meta

  if (nrow(labels_anc) == 0) {
    message("No PRS files for ancestry: ", code); next
  }
  if (nrow(meta_anc) == 0) {
    message("No individuals in meta (check META_PATH): ", code); next
  }

  message("\n=== Ancestry: ", anc_label,
          " | mz_mode=", MZ_MODE,
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
    apply_bh()

  out_path <- file.path(OUT_DIR, paste0("within_sibs_unaffected_", gsub(" ", "_", anc_label), ".tsv"))
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
      paste0("paper_ready_", paper_tag, "_sibs_unaffected_",
             gsub(" ", "_", anc_label), ".csv")
    )
    readr::write_csv(paper_anc, paper_path)
    message("Wrote paper-ready: ", paper_path)
  }

  all_anc_res[[anc_label]] <- res_anc
}

# =========================================
# EUR-only summary
# =========================================
res_all_summary <- bind_rows(all_anc_res)

cat("\n===== Unaffected-sibling EUR-only sample | MZ_MODE:", MZ_MODE,
    "| BH_SCOPE:", BH_SCOPE, "=====\n")

cat("\nSample (fitted models):\n")
print(as.data.frame(
  res_all_summary %>% filter(!is.na(estimate)) %>%
    summarise(n_ind = max(n_ind), n_fam = max(n_fam),
              models_fitted = n(),
              median_disc_fam = median(n_disc_fam, na.rm = TRUE))
))

cat("\nFDR-significant (q<0.05) by term:\n")
print(as.data.frame(
  res_all_summary %>% filter(!is.na(estimate)) %>%
    group_by(term) %>%
    summarise(fitted = n(),
              sig = sum(!is.na(q_BH) & q_BH < 0.05),
              median_abs_beta = median(abs(estimate), na.rm = TRUE),
              median_se = median(std.error, na.rm = TRUE),
              .groups = "drop")
))

cat("\nDropped models, by reason:\n")
print(as.data.frame(
  res_all_summary %>% filter(is.na(estimate)) %>% count(note, sort = TRUE)
))

fwrite(
  res_all_summary,
  file.path(OUT_DIR, "within_sibs_unaffected_EUR_all_results.tsv"),
  sep = "\t"
)

message("\nWrote EUR-only unaffected sibling results to: ", OUT_DIR)
