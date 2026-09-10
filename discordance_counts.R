# =========================================================
# SPARK V20 descriptive / discordance table
# Primary within-family sibling analytic sample
#
# =========================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
})

# ========= USER INPUTS =========
PRS_DIR        <- "/Users/asgelz01/Downloads/PRS_all/PRS"
ROLES_PATH     <- "/Users/asgelz01/Downloads/SPARKDataRelease_2026-06-25/roles-2026-06-25.csv"
INDIV_PATH     <- "/Users/asgelz01/Downloads/SPARKDataRelease_2026-06-25/individuals_registration-2026-06-25.csv"
BASIC_MED_PATH <- "/Users/asgelz01/Downloads/SPARKDataRelease_2026-06-25/basic_medical_screening-2026-06-25.csv"

IWES_META_PATH <- "/Users/asgelz01/Downloads/SPARK.iWES_v3.2024_08.sample_metadata.tsv"

OUT_PATH <- "/Users/asgelz01/Downloads/PRS_all/discordance_V20_primary_siblings.tsv"

outcomes <- c(
  "birth_etoh_subst",
  "birth_ivh",
  "birth_oxygen",
  "birth_pg_inf",
  "birth_prem",
  "growth_low_wt",
  "growth_macroceph",
  "growth_microceph",
  "med_cond_birth",
  "med_cond_birth_def",
  "med_cond_growth"
)

TARGET_ANCES <- c("AFR", "AMR", "EUR")

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

# =========================================================
# 1. Ancestry map from PRS files
# =========================================================
sscore_files <- list.files(PRS_DIR, pattern = "\\.sscore$", full.names = TRUE)
if (length(sscore_files) == 0) stop("No .sscore files found in: ", PRS_DIR)

label_from_sscore <- function(fp) {
  b <- basename(fp)

  m1 <- stringr::str_match(b, "^([A-Za-z0-9]+)_(AFR|AMR)\\.sscore$")
  if (!is.na(m1[1,1])) {
    return(tibble(trait = m1[1,2], ancestry = m1[1,3], file = fp))
  }

  m2 <- stringr::str_match(b, "^EUR_([A-Za-z0-9]+)\\.sscore$")
  if (!is.na(m2[1,1])) {
    return(tibble(trait = m2[1,2], ancestry = "EUR", file = fp))
  }

  NULL
}

labels <- map_dfr(sscore_files, label_from_sscore) %>%
  filter(ancestry %in% TARGET_ANCES)

if (nrow(labels) == 0) stop("No ancestry-specific PRS files matched expected names.")

read_sscore_iids <- function(path) {
  dt <- fread(path, sep = "\t", header = TRUE, check.names = FALSE)
  setnames(dt, sub("^#", "", names(dt)))
  if (!"IID" %in% names(dt)) stop("IID missing in ", path)
  tibble(IID = as.character(dt$IID))
}

anc_map <- labels %>%
  mutate(ids = map(file, read_sscore_iids)) %>%
  unnest(ids) %>%
  select(IID, ancestry) %>%
  distinct()

dup_anc <- anc_map %>% count(IID) %>% filter(n > 1)
if (nrow(dup_anc) > 0) {
  warning(nrow(dup_anc), " IIDs occur in >1 ancestry PRS set; dropping them.")
  anc_map <- anc_map %>% filter(!IID %in% dup_anc$IID)
}

# =========================================================
# 2. V20 explicit proband/sibling sample definition
# =========================================================
blank_na <- function(x) {
  x <- trimws(as.character(x))
  if_else(x %in% c("", "NA", "0", "."), NA_character_, x)
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

role_subject_ids     <- unique(roles$IID)
affected_sibling_ids <- split_ids(roles$affected_sibs)
control_sibling_ids  <- split_ids(roles$control_sibs)

child_ids <- unique(c(
  role_subject_ids,
  affected_sibling_ids,
  control_sibling_ids
))

parent_ids <- unique(c(roles$biomother, roles$biofather))
parent_ids <- parent_ids[!is.na(parent_ids)]
parent_child_overlap <- intersect(child_ids, parent_ids)

indiv <- fread(INDIV_PATH, sep = ",", header = TRUE, quote = "\"") %>%
  as_tibble() %>%
  transmute(
    IID = as.character(subject_sp_id),
    FID = as.character(family_sf_id)
  ) %>%
  distinct(IID, .keep_all = TRUE)

sibling_ids <- unique(c(affected_sibling_ids, control_sibling_ids))
missing_sibling_ids <- setdiff(sibling_ids, indiv$IID)

message("V20 explicit proband/sibling construction:")
message("  roles-row IDs: ", length(role_subject_ids))
message("  affected sibling IDs: ", length(affected_sibling_ids))
message("  control sibling IDs: ", length(control_sibling_ids))
message("  unique child IDs total: ", length(child_ids))
message("  sibling IDs absent from individuals_registration: ",
        length(missing_sibling_ids))
message("  explicit child IDs also listed as biological parents: ",
        length(parent_child_overlap))

if (length(missing_sibling_ids) > 0) {
  warning(
    "Some sibling IDs listed in roles are absent from individuals_registration. Example(s): ",
    paste(head(missing_sibling_ids, 10), collapse = ", ")
  )
}

child_meta <- indiv %>%
  filter(IID %in% child_ids, !is.na(FID)) %>%
  distinct(IID, .keep_all = TRUE)

# =========================================================
# 3. MZ sets: same SPID-linked logic as primary within-family model
# =========================================================
iwes <- fread(IWES_META_PATH, sep = "\t", header = TRUE) %>%
  as_tibble() %>%
  transmute(
    IID = as.character(spid),
    identical_twins = as.character(identical_twins)
  )

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

  if (!nrow(edges))
    return(tibble(IID = character(0), mz_group = character(0)))

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
    ra <- find_root(edges$a[k])
    rb <- find_root(edges$b[k])
    if (ra != rb) parent[[rb]] <- ra
  }

  tibble(
    IID = nodes,
    mz_group = vapply(nodes, find_root, character(1), USE.NAMES = FALSE)
  )
}

mz_groups <- build_mz_groups(iwes$IID, iwes$identical_twins)

keep_one_per_mz_group <- function(df) {
  if (!nrow(df)) return(df)

  df %>%
    group_by(FID, mz_group) %>%
    filter(is.na(mz_group) | IID == min(IID)) %>%
    ungroup()
}

message("MZ linkage: ", nrow(mz_groups), " linked individuals in ",
        n_distinct(mz_groups$mz_group), " MZ sets.")

# =========================================================
# 4. V20 Basic Medical Screening coding
# =========================================================
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

basic_med <- fread(
  BASIC_MED_PATH,
  sep = ",",
  header = TRUE,
  quote = "\"",
  colClasses = "character"
) %>%
  as_tibble() %>%
  select(subject_sp_id, all_of(outcomes)) %>%
  rename(IID = subject_sp_id) %>%
  mutate(
    IID = as.character(IID),
    across(all_of(outcomes), clean_v20_binary)
  ) %>%
  distinct(IID, .keep_all = TRUE)

# =========================================================
# 5. Primary sibling sample
# =========================================================
df <- child_meta %>%
  inner_join(anc_map, by = "IID") %>%
  inner_join(basic_med, by = "IID") %>%
  left_join(mz_groups, by = "IID")

message("Before MZ dedup: ", nrow(df), " genotyped children in ",
        n_distinct(df$FID), " V20 families.")

n_before_mz <- nrow(df)
df <- keep_one_per_mz_group(df)
message("MZ deduplication removed ", n_before_mz - nrow(df),
        " individual(s); one SPID retained per MZ set.")

# =========================================================
# 6. Outcome-specific individual descriptives
# =========================================================
descriptives_one <- function(dat, cond) {
  dat %>%
    filter(!is.na(.data[[cond]])) %>%
    group_by(ancestry) %>%
    summarise(
      individuals_with_data = n(),
      families_with_data    = n_distinct(FID),
      events                 = sum(.data[[cond]] == 1L),
      non_events             = sum(.data[[cond]] == 0L),
      .groups = "drop"
    ) %>%
    mutate(condition = cond) %>%
    select(ancestry, condition, individuals_with_data,
           families_with_data, events, non_events)
}

desc_tbl <- map_dfr(outcomes, ~ descriptives_one(df, .x))

# =========================================================
# 7. Outcome-specific within-family discordance
# =========================================================
summarize_family_cond <- function(dat, cond) {
  dat %>%
    filter(!is.na(.data[[cond]])) %>%
    group_by(ancestry, FID) %>%
    summarise(
      n  = n(),
      n1 = sum(.data[[cond]] == 1L),
      n0 = sum(.data[[cond]] == 0L),
      .groups = "drop"
    ) %>%
    mutate(
      eligible_family   = as.integer(n >= 2),
      total_pairs       = ifelse(n >= 2, choose(n, 2), 0),
      discordant_pairs  = ifelse(n >= 2, n0 * n1, 0),
      discordant_family = as.integer(n >= 2 & n1 > 0 & n0 > 0)
    )
}

disc_tbl <- map_dfr(outcomes, function(cond) {
  fam <- summarize_family_cond(df, cond)

  fam %>%
    group_by(ancestry) %>%
    summarise(
      discordant_families = sum(discordant_family),
      total_families      = sum(eligible_family),
      discordant_siblings = sum(discordant_pairs),
      total_sibling_pairs = sum(total_pairs),
      events_in_disc_fams = sum(ifelse(discordant_family == 1L, n1, 0)),
      .groups = "drop"
    ) %>%
    mutate(condition = cond) %>%
    select(
      ancestry, condition,
      discordant_families, total_families,
      discordant_siblings, total_sibling_pairs,
      events_in_disc_fams
    )
})

# =========================================================
# 8. Final table
# =========================================================
final_tbl <- desc_tbl %>%
  left_join(disc_tbl, by = c("ancestry", "condition")) %>%
  mutate(
    Condition = recode(condition, !!!outcome_map, .default = condition)
  ) %>%
  arrange(ancestry, Condition) %>%
  transmute(
    Ancestry = ancestry,
    Condition,
    `Individuals with outcome data` = coalesce(as.integer(individuals_with_data), 0L),
    `Affected individuals`          = coalesce(as.integer(events), 0L),
    `Unaffected individuals`        = coalesce(as.integer(non_events), 0L),
    `Families (>=2 siblings with outcome data)` =
      coalesce(as.integer(total_families), 0L),
    `Discordant sibling pairs`      = coalesce(as.integer(discordant_siblings), 0L),
    `Discordant families`           = coalesce(as.integer(discordant_families), 0L),
    `Events in discordant families` =
      coalesce(as.integer(events_in_disc_fams), 0L)
  )

fwrite(final_tbl, OUT_PATH, sep = "\t", quote = FALSE)
message("Wrote: ", OUT_PATH)
print(final_tbl, n = nrow(final_tbl))
