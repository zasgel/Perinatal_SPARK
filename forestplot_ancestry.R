# =========================================================
# Forest plots for population-level results by ancestry
#
# Compares:
#   - European
#   - African
#   - Admixed American
#
# =========================================================

# ---------- 1) Packages ----------
packs <- c("tidyverse", "scales", "patchwork")

to_install <- setdiff(packs, rownames(installed.packages()))
if (length(to_install) > 0) {
  install.packages(to_install, repos = "https://cloud.r-project.org")
}
invisible(lapply(packs, library, character.only = TRUE))

# ---------- 2) Paths ----------
base_dir <- "/Users/asgelz01/Downloads/PRS_all/PRS_ASD/per_ancestry_outputs_new"

eur_path <- file.path(base_dir, "prs_perinatal_European.csv")
afr_path <- file.path(base_dir, "prs_perinatal_African.csv")
amr_path <- file.path(base_dir, "prs_perinatal_Admixed_American.csv")

# ---------- 3) Outputs ----------
out_dir_ind   <- "/Users/asgelz01/Downloads/PRS_all/plots_population_ancestry_individual"
out_dir_panel <- "/Users/asgelz01/Downloads/PRS_all/plots_population_ancestry"

dir.create(out_dir_ind, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_panel, showWarnings = FALSE, recursive = TRUE)

# ---------- 4) Trait/outcome settings ----------
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

preferred_outcome_order <- c(
  "Birth defects",
  "Birth or pregnancy complications",
  "Difficulty gaining weight",
  "Fetal alcohol syndrome or in-utero drug/alcohol exposure",
  "Growth abnormalities",
  "Intraventricular hemorrhage",
  "Macrocephaly",
  "Microcephaly",
  "Oxygen deprivation at birth requiring NICU stay",
  "Premature birth",
  "Serious prenatal infection"
)

# ---------- 5) Helpers ----------
empty_long <- tibble(
  Trait    = character(),
  Outcome  = character(),
  Ancestry = character(),
  OR       = numeric(),
  CI_low   = numeric(),
  CI_high  = numeric(),
  p        = numeric()
)

safe_read_csv <- function(path, ...) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) {
    message("Missing file: ", path)
    return(NULL)
  }
  
  tryCatch(
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, ...),
    error = function(e) {
      message("Could not read: ", path)
      message("Reason: ", e$message)
      NULL
    }
  )
}

safe_recode <- function(x, map_vec) {
  dplyr::recode(x, !!!map_vec, .default = x)
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

sanitize_filename <- function(x) {
  x |>
    gsub("[^A-Za-z0-9]+", "_", x = _) |>
    gsub("^_|_$", "", x = _)
}

# ---------- 6) Reader / harmonizer ----------
coerce_population_tbl <- function(tbl, ancestry_name) {
  if (is.null(tbl)) return(empty_long)
  
  has_or_col   <- "Odds ratio" %in% names(tbl)
  has_beta_col <- "Beta" %in% names(tbl)
  
  if (!"Trait" %in% names(tbl) || !"Outcome" %in% names(tbl)) {
    message("File missing Trait or Outcome column.")
    return(empty_long)
  }
  
  if (has_or_col) {
    out <- tbl %>%
      transmute(
        Trait    = harmonize_trait_names(as.character(Trait)),
        Outcome  = safe_recode(as.character(Outcome), outcome_map),
        Ancestry = ancestry_name,
        OR       = as.numeric(`Odds ratio`),
        CI_low   = as.numeric(`CI lower`),
        CI_high  = as.numeric(`CI upper`),
        p        = if ("p-value" %in% names(.)) as.numeric(`p-value`) else NA_real_
      )
  } else if (has_beta_col) {
    out <- tbl %>%
      transmute(
        Trait    = harmonize_trait_names(as.character(Trait)),
        Outcome  = safe_recode(as.character(Outcome), outcome_map),
        Ancestry = ancestry_name,
        OR       = exp(as.numeric(Beta)),
        CI_low   = as.numeric(`CI lower`),
        CI_high  = as.numeric(`CI upper`),
        p        = if ("p-value" %in% names(.)) as.numeric(`p-value`) else NA_real_
      )
  } else {
    message("Population file missing expected Odds ratio or Beta column.")
    return(empty_long)
  }
  
  out
}

# ---------- 7) Read files ----------
eur_raw <- safe_read_csv(eur_path)
afr_raw <- safe_read_csv(afr_path)
amr_raw <- safe_read_csv(amr_path)

# ---------- 8) Convert to plotting table ----------
eur <- coerce_population_tbl(eur_raw, "European")
afr <- coerce_population_tbl(afr_raw, "African")
amr <- coerce_population_tbl(amr_raw, "Admixed American")

dat <- bind_rows(eur, afr, amr) %>%
  filter(
    is.finite(OR),
    is.finite(CI_low),
    is.finite(CI_high),
    !is.na(Trait),
    !is.na(Outcome),
    !is.na(Ancestry)
  ) %>%
  filter(
    !Trait %in% c(
      "Gestational diabetes",
      "Postpartum hemorrhage",
      "Pre-eclampsia"
    )
  ) %>%
  mutate(
    Ancestry = factor(
      Ancestry,
      levels = c("European", "African", "Admixed American")
    ),
    Trait_disp = Trait,
    Outcome_disp = Outcome
  )

if (nrow(dat) == 0) {
  stop("No valid rows available for plotting after harmonization.")
}

# ---------- 9) Enforce outcome order ----------
present_outcomes <- preferred_outcome_order[preferred_outcome_order %in% unique(dat$Outcome_disp)]
other_outcomes   <- setdiff(sort(unique(dat$Outcome_disp)), preferred_outcome_order)
outcome_levels   <- rev(c(present_outcomes, other_outcomes))

dat <- dat %>%
  mutate(
    Outcome_disp = factor(Outcome_disp, levels = outcome_levels)
  )

# ---------- 10) Plot function ----------
make_trait_plot <- function(df_trait, title_str, x_pad = 0.10, x_min = 0.5, x_max = 3) {
  
  rng <- range(c(df_trait$CI_low, df_trait$CI_high), na.rm = TRUE)
  
  lo <- max(x_min, 10^(log10(rng[1]) - x_pad))
  hi <- min(x_max, 10^(log10(rng[2]) + x_pad))
  
  df_trait <- df_trait %>%
    mutate(
      CI_low_capped  = pmax(CI_low, lo),
      CI_high_capped = pmin(CI_high, hi)
    )
  
  pd <- position_dodge(width = 0.60)
  
  ggplot(df_trait, aes(x = OR, y = Outcome_disp, color = Ancestry, shape = Ancestry)) +
    geom_vline(xintercept = 1, linetype = 2) +
    geom_errorbarh(
      aes(xmin = CI_low_capped, xmax = CI_high_capped),
      position = pd,
      height = 0.25,
      linewidth = 0.4
    ) +
    geom_point(position = pd, size = 2.3) +
    scale_x_log10(
      limits = c(lo, hi),
      breaks = c(0.5, 0.67, 1, 1.5, 2, 3),
      labels = scales::label_number(accuracy = 0.01)
    ) +
    scale_shape_manual(
      values = c(
        "European"         = 16,
        "African"          = 17,
        "Admixed American" = 4
      )
    ) +
    labs(
      title = title_str,
      x = "Odds ratio (log scale)",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      axis.text.y = element_text(size = 11),
      legend.title = element_blank()
    )
}

# ---------- 11) Build individual plots ----------
traits <- sort(unique(dat$Trait_disp))

plots <- setNames(vector("list", length(traits)), traits)

for (tr in traits) {
  df_tr <- dat %>% filter(Trait_disp == tr)
  
  p <- make_trait_plot(df_tr, tr)
  plots[[tr]] <- p
  
  safe_name <- sanitize_filename(tr)
  
  ggsave(
    filename = file.path(out_dir_ind, paste0("forest_population_ancestry_", safe_name, ".png")),
    plot = p,
    width = 6,
    height = 6.5,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(out_dir_ind, paste0("forest_population_ancestry_", safe_name, ".pdf")),
    plot = p,
    width = 6,
    height = 6.5
  )
}

# ---------- 12) Panel pages ----------
chunk_vec <- function(x, n) {
  split(x, ceiling(seq_along(x) / n))
}

page_chunks <- chunk_vec(traits, 4)
panel_list <- list()

for (i in seq_along(page_chunks)) {
  tr_page <- page_chunks[[i]]
  
  page_plots <- lapply(seq_along(tr_page), function(j) {
    tr <- tr_page[j]
    p  <- plots[[tr]]
    
    if (j %% 2 == 0) {
      p <- p + theme(
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank()
      )
    }
    p
  })
  
  if (length(page_plots) < 4) {
    page_plots <- c(
      page_plots,
      replicate(4 - length(page_plots), ggplot() + theme_void(), simplify = FALSE)
    )
  }
  
  panel_i <- wrap_plots(page_plots, ncol = 2) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  ggsave(
    filename = file.path(out_dir_panel, sprintf("forest_population_ancestry_page_%02d.png", i)),
    plot = panel_i,
    width = 14,
    height = 10,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(out_dir_panel, sprintf("forest_population_ancestry_page_%02d.pdf", i)),
    plot = panel_i,
    width = 14,
    height = 10
  )
  
  panel_list[[i]] <- panel_i
}

# ---------- 13) Combined PDF ----------
pdf_path <- file.path(out_dir_panel, "forest_population_ancestry_ALL.pdf")

pdf(pdf_path, width = 14, height = 20)
for (p in panel_list) {
  print(p)
}
dev.off()

message("Done.")
message("Individual plots saved to: ", out_dir_ind)
message("Panel plots saved to: ", out_dir_panel)
message("Combined PDF saved to: ", pdf_path)

