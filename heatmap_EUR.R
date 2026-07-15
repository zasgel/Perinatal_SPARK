# ================================
# Heatmap of European-only PRS Effects
# Population-level + Within-family
# ================================

packs <- c("tidyverse", "scales")
invisible(lapply(packs, library, character.only = TRUE))

# ---------- Paths ----------
eur_population_path <- "/Users/asgelz01/Downloads/PRS_all/PRS_ASD/per_ancestry_outputs_new/prs_perinatal_European.csv"
eur_within_path <- "/Users/asgelz01/Downloads/PRS_all/within_sibs_results_per_ancestry/paper_ready_within_sibs_European.csv"

out_dir_root    <- "/Users/asgelz01/Downloads/PRS_all/plots_heatmaps_EUR"
out_dir_figures <- file.path(out_dir_root, "figures")

dir.create(out_dir_figures, showWarnings = FALSE, recursive = TRUE)

# ---------- Outcome order ----------
preferred_outcome_order <- c(
  "Birth defects",
  "Birth or pregnancy complications",
  "Difficulty gaining weight",
  "Fetal alcohol or substance exposure",
  "Growth abnormalities",
  "Intraventricular hemorrhage",
  "Macrocephaly",
  "Microcephaly",
  "Oxygen deprivation at birth requiring NICU stay",
  "Premature birth",
  "Serious prenatal infection"
)

# ---------- Helpers ----------
empty_long <- tibble(
  Trait   = character(),
  Outcome = character(),
  Model   = character(),
  OR      = numeric(),
  CI_low  = numeric(),
  CI_high = numeric(),
  p       = numeric(),
  q       = numeric()
)

safe_read_csv <- function(path, ...) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) {
    message("Missing file: ", path)
    return(NULL)
  }
  
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, ...)
}

harmonize_trait_names <- function(x) {
  case_when(
    x %in% c("AS", "Asthma") ~ "Asthma",
    x %in% c("OB", "Obesity") ~ "Obesity",
    x %in% c("MDD", "Major depressive disorder", "Major depression") ~ "Major depression",
    x %in% c("SCZ", "Schizophrenia") ~ "Schizophrenia",
    x %in% c("PE", "Pre-eclampsia", "Preeclampsia") ~ "Pre-eclampsia",
    x %in% c("PPH", "Postpartum hemorrhage") ~ "Postpartum hemorrhage",
    x %in% c("GD", "Gestational diabetes") ~ "Gestational diabetes",
    TRUE ~ x
  )
}

coerce_population_tbl <- function(tbl, model_name) {
  if (is.null(tbl)) return(empty_long)
  
  tbl %>%
    transmute(
      Trait   = harmonize_trait_names(as.character(Trait)),
      Outcome = as.character(Outcome),
      Model   = model_name,
      OR      = exp(as.numeric(Beta)),
      CI_low  = as.numeric(`CI lower`),
      CI_high = as.numeric(`CI upper`),
      p       = as.numeric(`p-value`),
      q       = as.numeric(`BH q-value`)
    )
}

coerce_within_tbl <- function(tbl, model_name) {
  if (is.null(tbl)) return(empty_long)
  
  tbl %>%
    transmute(
      Trait   = harmonize_trait_names(as.character(Trait)),
      Outcome = as.character(Outcome),
      Model   = model_name,
      OR      = as.numeric(`Odds ratio`),
      CI_low  = as.numeric(`CI lower`),
      CI_high = as.numeric(`CI upper`),
      p       = as.numeric(`p-value`),
      q       = as.numeric(`BH q-value`)
    )
}

# ---------- Read files ----------
eur_population_raw <- safe_read_csv(eur_population_path)
eur_within_raw     <- safe_read_csv(eur_within_path)

# ---------- Convert ----------
population <- coerce_population_tbl(
  eur_population_raw,
  "Population-level"
)

withinfam <- coerce_within_tbl(
  eur_within_raw,
  "Within-family"
)

# ---------- Combine ----------
dat <- bind_rows(
  population,
  withinfam
) %>%
  mutate(
    Outcome = recode(
      Outcome,
      "Fetal alcohol syndrome or in-utero drug/alcohol exposure" =
        "Fetal alcohol or substance exposure"
    )
  ) %>%
  filter(
    Outcome != "Serious prenatal infection"
  ) %>%
  filter(
    is.finite(OR),
    OR > 0,
    is.finite(CI_low),
    is.finite(CI_high),
    !is.na(Trait),
    !is.na(Outcome),
    !is.na(Model)
  ) %>%
  mutate(
    Trait = harmonize_trait_names(Trait),
    Model = factor(Model, levels = c("Population-level", "Within-family"))
  )

if (nrow(dat) == 0) {
  stop("No valid rows available for plotting.")
}

# ---------- Ordering ----------
trait_levels <- sort(unique(dat$Trait))

present_outcomes <- preferred_outcome_order[
  preferred_outcome_order %in% unique(dat$Outcome)
]

other_outcomes <- setdiff(
  sort(unique(dat$Outcome)),
  preferred_outcome_order
)

outcome_levels <- rev(c(present_outcomes, other_outcomes))

dat <- dat %>%
  mutate(
    Trait = factor(Trait, levels = trait_levels),
    Outcome = factor(Outcome, levels = outcome_levels)
  )

# ---------- Direction consistency ----------
agree_df <- dat %>%
  select(Trait, Outcome, Model, OR) %>%
  pivot_wider(names_from = Model, values_from = OR) %>%
  filter(
    is.finite(`Population-level`),
    `Population-level` > 0,
    is.finite(`Within-family`),
    `Within-family` > 0
  ) %>%
  mutate(
    sign_population = sign(log(`Population-level`)),
    sign_within     = sign(log(`Within-family`)),
    direction_consistent = sign_population == sign_within & sign_population != 0
  )

# ---------- Heatmap data ----------
hm <- dat %>%
  filter(Model == "Population-level") %>%
  left_join(
    agree_df %>%
      select(Trait, Outcome, direction_consistent),
    by = c("Trait", "Outcome")
  ) %>%
  mutate(
    sig_label = case_when(
      !is.na(q) & q < 0.05 ~ "**",
      !is.na(p) & p < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# ---------- Theme ----------
theme_hm <- theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.title = element_text(size = 11)
  )

# ---------- Heatmap 1: significance stars ----------
p_sig <- ggplot(hm, aes(x = Trait, y = Outcome, fill = log(OR))) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0,
    name = "log(OR)"
  ) +
  geom_text(
    aes(label = sig_label),
    size = 5,
    color = "black",
    fontface = "bold"
  ) +
  labs(
    title = "European ancestry population-level PRS associations with perinatal outcomes",
    subtitle = "Fill = log(OR); * p < 0.05; ** BH q < 0.05",
    x = "Trait",
    y = "Outcome"
  ) +
  theme_hm

# ---------- Heatmap 2: direction consistency ----------
p_direction <- ggplot(hm, aes(x = Trait, y = Outcome, fill = log(OR))) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0,
    name = "log(OR)"
  ) +
  geom_text(
    data = hm %>% filter(direction_consistent %in% TRUE),
    label = "*",
    size = 6,
    color = "black",
    fontface = "bold"
  ) +
  labs(
    title = "European ancestry population-level PRS associations with perinatal outcomes",
    subtitle = "Fill = log(OR); * indicates concordant effect direction in European within-family analyses",
    x = "Trait",
    y = "Outcome"
  ) +
  theme_hm

# ---------- Save ----------
ggsave(
  file.path(out_dir_figures, "heatmap_EUR_population_significance.png"),
  p_sig,
  width = 12,
  height = 8.5,
  dpi = 300
)

ggsave(
  file.path(out_dir_figures, "heatmap_EUR_population_significance.pdf"),
  p_sig,
  width = 12,
  height = 8.5
)

ggsave(
  file.path(out_dir_figures, "heatmap_EUR_population_direction_consistency.png"),
  p_direction,
  width = 12,
  height = 8.5,
  dpi = 300
)

ggsave(
  file.path(out_dir_figures, "heatmap_EUR_population_direction_consistency.pdf"),
  p_direction,
  width = 12,
  height = 8.5
)

# ---------- Export consistency table ----------
consistent_csv <- file.path(out_dir_root, "EUR_direction_consistent_pairs.csv")

agree_df %>%
  mutate(
    consistent_flag = if_else(direction_consistent, "consistent", "not_consistent")
  ) %>%
  arrange(Outcome, Trait) %>%
  write.csv(consistent_csv, row.names = FALSE)

message("Done.")
message("Saved EUR-only heatmaps to: ", out_dir_figures)
message("Saved consistency table to: ", consistent_csv)
