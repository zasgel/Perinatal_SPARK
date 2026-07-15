
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
})

# ========= USER INPUTS =========
PC_DIR          <- "/Users/asgelz01/Downloads/PRS_all/PRS"
META_PATH       <- "/Users/asgelz01/Downloads/SPARK.iWES_v3.2024_08.sample_metadata.tsv"
BASIC_MED_PATH  <- "/Users/asgelz01/Downloads/SPARKDataRelease_2025-07-14/basic_medical_screening-2025-07-14.csv"
OUT_PATH        <- "discordance_AMR_EUR_AFR_with_individual_counts.tsv"

# Outcomes / conditions to include
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

TARGET_ANCES <- c("AFR", "AMR", "EUR")   # restrict to these ancestries
EXCLUDE_IDENTICAL_TWINS <- FALSE         # set TRUE to drop identical twins

# ========= BUILD ANCESTRY MAP FROM THE SAME PRS INPUTS  =========

sscore_files <- list.files(PC_DIR, pattern = "\\.sscore$", full.names = TRUE)
if (length(sscore_files) == 0) stop("No .sscore files found in: ", PC_DIR)

label_from_sscore <- function(fp) {
  b <- basename(fp)
  
  # AFR/AMR style: AS_AFR.sscore, MDD_AMR.sscore, OB_AFR.sscore
  m1 <- stringr::str_match(b, "^([A-Z]+)_(AFR|AMR)\\.sscore$")
  if (!is.na(m1[1,1])) {
    return(list(method = "nonEUR", trait = m1[1,2], ancestry = m1[1,3], file = fp))
  }
  
  m2 <- stringr::str_match(b, "^EUR_([A-Z]+)\\.sscore$")
  if (!is.na(m2[1,1])) {
    return(list(method = "EUR", trait = m2[1,2], ancestry = "EUR", file = fp))
  }
  
  # ignore META files and anything else
  return(list(method = "Unknown", trait = b, ancestry = NA_character_, file = fp))
}

labels <- purrr::map(sscore_files, label_from_sscore) %>% dplyr::bind_rows()
labels <- labels %>% dplyr::filter(!is.na(ancestry), ancestry %in% TARGET_ANCES)

read_sscore_iids <- function(path) {
  dt <- data.table::fread(path, sep = "\t", header = TRUE, check.names = FALSE)
  data.table::setnames(dt, sub("^#", "", names(dt)))
  if (!"IID" %in% names(dt)) stop("IID missing in ", path)
  tibble::tibble(IID = as.character(dt$IID))
}

anc_map <- labels %>%
  dplyr::mutate(ids = purrr::map(file, read_sscore_iids)) %>%
  tidyr::unnest(ids) %>%
  dplyr::select(IID, ancestry) %>%
  dplyr::distinct()



meta <- data.table::fread(META_PATH, sep = "\t", quote = "", header = TRUE) |> 
  tibble::as_tibble()


meta <- meta |>
  mutate(
    identical_twins_raw = trimws(as.character(identical_twins)),
    twin_partner = dplyr::if_else(stringr::str_detect(identical_twins_raw, "^SP\\d+$"),
                                  identical_twins_raw, NA_character_),
    identical_twins = case_when(
      !is.na(twin_partner) ~ TRUE,
      tolower(identical_twins_raw) %in% c("true","1","yes","y") ~ TRUE,
      TRUE ~ FALSE
    )
  )

basic_med <- fread(BASIC_MED_PATH) %>%
  select(subject_sp_id, dplyr::all_of(outcomes)) %>%
  rename(IID = subject_sp_id) %>%
  mutate(across(all_of(outcomes), ~ tidyr::replace_na(., 0L))) %>%
  mutate(across(all_of(outcomes), as.integer))

# ========= JOIN: IID -> ancestry (PCs) and family (meta) =========
df <- basic_med %>%
  dplyr::left_join(anc_map, by = c("IID" = "IID")) %>%
  dplyr::rename(spid = IID) %>%
  dplyr::left_join(meta, by = "spid") %>%
  dplyr::mutate(
    FID = sfid,
    ancestry = ancestry
  ) %>%
  dplyr::select(spid, FID, ancestry, identical_twins, twin_partner, dplyr::all_of(outcomes)) %>%
  dplyr::filter(!is.na(FID), !is.na(ancestry), ancestry %in% TARGET_ANCES)

if (EXCLUDE_IDENTICAL_TWINS) {
  df <- df %>% dplyr::filter(!identical_twins)
}

keep_one_per_mz_pair <- function(dd) {
  if (!all(c("spid","FID","twin_partner") %in% names(dd))) return(dd)
  dd |>
    mutate(tp_key = dplyr::case_when(
      is.na(twin_partner) ~ NA_character_,
      spid < twin_partner ~ paste0(spid, "__", twin_partner),
      TRUE                ~ paste0(twin_partner, "__", spid)
    )) |>
    arrange(FID, tp_key, spid) |>
    group_by(tp_key) |>
    # keep one per MZ pair; keep all non-twins
    filter(is.na(tp_key) | dplyr::row_number() == 1L) |>
    ungroup() |>
    select(-tp_key)
}

df <- df %>% keep_one_per_mz_pair()


# ========= INDIVIDUAL-LEVEL DESCRIPTIVES =========
descriptives_one <- function(dat, cond) {
  dat %>%
    group_by(ancestry) %>%
    summarise(
      individuals_n = n(),
      families_n    = n_distinct(FID),
      events        = sum(.data[[cond]] == 1L),
      non_events    = sum(.data[[cond]] == 0L),
      .groups = "drop"
    ) %>%
    mutate(condition = cond) %>%
    select(ancestry, condition, individuals_n, families_n, events, non_events)
}

desc_list <- lapply(outcomes, function(cn) descriptives_one(df, cn))
desc_tbl  <- bind_rows(desc_list)

# ========= WITHIN-FAMILY DISCORDANCE =========
summarize_family_cond <- function(dat, cond) {
  dat %>%
    group_by(ancestry, FID) %>%
    summarise(
      n  = dplyr::n(),
      n1 = sum(.data[[cond]] == 1L),
      n0 = sum(.data[[cond]] == 0L),
      min_val = min(.data[[cond]]),
      max_val = max(.data[[cond]]),
      .groups = "drop"
    ) %>%
    mutate(
      eligible_family   = as.integer(n >= 2),
      total_pairs       = ifelse(n >= 2, choose(n, 2), 0),
      discordant_pairs  = ifelse(n >= 2, n0 * n1, 0),
      discordant_family = ifelse(n >= 2, as.integer(min_val != max_val), 0)
    )
}

discord_list <- lapply(outcomes, function(cn) {
  fam_summ <- summarize_family_cond(df, cn)
  fam_summ %>%
    group_by(ancestry) %>%
    summarise(
      discordant_families = sum(discordant_family),
      total_families      = sum(eligible_family),
      discordant_siblings = sum(discordant_pairs),
      total_sibling_pairs = sum(total_pairs),
      .groups = "drop"
    ) %>%
    mutate(condition = cn) %>%
    select(ancestry, condition,
           discordant_families, total_families,
           discordant_siblings, total_sibling_pairs)
})
disc_tbl <- bind_rows(discord_list)


final_tbl <- desc_tbl %>%
  dplyr::left_join(disc_tbl, by = c("ancestry","condition")) %>%
  dplyr::arrange(ancestry, condition)

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

final_tbl <- final_tbl %>%
  dplyr::mutate(
    condition = dplyr::recode(condition, !!!outcome_map, .default = condition)
  ) %>%
  dplyr::transmute(
    Ancestry = ancestry,
    Condition = condition,
    `Families (≥2 siblings)` = coalesce(as.integer(total_families), 0L),
    `Affected individuals` = coalesce(as.integer(events), 0L),
    `Discordant sibling pairs` = coalesce(as.integer(discordant_siblings), 0L)
  )

data.table::fwrite(final_tbl, OUT_PATH, sep = "\t", quote = FALSE)
message("Wrote: ", OUT_PATH)
print(final_tbl, n = nrow(final_tbl))
