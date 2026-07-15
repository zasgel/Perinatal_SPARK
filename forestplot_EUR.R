# =========================================================
# EUR-only forest plots for:
# Population-level, Between-family, Within-family
# =========================================================

packs <- c("tidyverse", "scales", "patchwork")
invisible(lapply(packs, library, character.only = TRUE))

# ---------- Paths ----------
eur_population_path <- "/Users/asgelz01/Downloads/PRS_all/PRS_ASD/per_ancestry_outputs_new/prs_perinatal_European.csv"

eur_within_path  <- "/Users/asgelz01/Downloads/PRS_all/within_sibs_results_per_ancestry/paper_ready_within_sibs_European.csv"
eur_between_path <- "/Users/asgelz01/Downloads/PRS_all/within_sibs_results_per_ancestry/paper_ready_between_sibs_European.csv"

out_dir_ind   <- "/Users/asgelz01/Downloads/PRS_all/plots_individual_EUR"
out_dir_panel <- "/Users/asgelz01/Downloads/PRS_all/plots_EUR"

dir.create(out_dir_ind, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_panel, showWarnings = FALSE, recursive = TRUE)

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
  "Premature birth"
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
  
  tryCatch(
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, ...),
    error = function(e) {
      message("Could not read: ", path)
      message("Reason: ", e$message)
      NULL
    }
  )
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

shorten_outcomes <- function(x) {
  dplyr::recode(
    x,
    "Fetal alcohol syndrome or in-utero drug/alcohol exposure" =
      "Fetal alcohol or substance exposure",
    .default = x
  )
}

sanitize_filename <- function(x) {
  x |>
    gsub("[^A-Za-z0-9]+", "_", x = _) |>
    gsub("^_|_$", "", x = _)
}

# ---------- Coercion functions ----------
coerce_population_tbl <- function(tbl, model_name) {
  if (is.null(tbl)) return(empty_long)
  
  tbl %>%
    transmute(
      Trait   = harmonize_trait_names(as.character(Trait)),
      Outcome = shorten_outcomes(as.character(Outcome)),
      Model   = model_name,
      OR      = if ("Odds ratio" %in% names(.)) {
        as.numeric(`Odds ratio`)
      } else {
        exp(as.numeric(Beta))
      },
      CI_low  = as.numeric(`CI lower`),
      CI_high = as.numeric(`CI upper`),
      p       = as.numeric(`p-value`),
      q       = as.numeric(`BH q-value`)
    )
}

coerce_sibling_tbl <- function(tbl, model_name) {
  if (is.null(tbl)) return(empty_long)
  
  tbl %>%
    transmute(
      Trait   = harmonize_trait_names(as.character(Trait)),
      Outcome = shorten_outcomes(as.character(Outcome)),
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
eur_between_raw    <- safe_read_csv(eur_between_path)

# ---------- Convert ----------
population <- coerce_population_tbl(
  eur_population_raw,
  "Population-level"
)

betweenfam <- coerce_sibling_tbl(
  eur_between_raw,
  "Between-family"
)

withinfam <- coerce_sibling_tbl(
  eur_within_raw,
  "Within-family"
)

# ---------- Combine ----------
dat <- bind_rows(
  population,
  betweenfam,
  withinfam
) %>%
  filter(
    Outcome != "Serious prenatal infection"
  ) %>%
  filter(
    is.finite(OR),
    is.finite(CI_low),
    is.finite(CI_high),
    OR > 0,
    CI_low > 0,
    CI_high > 0,
    !is.na(Trait),
    !is.na(Outcome),
    !is.na(Model)
  ) %>%
  mutate(
    Trait = harmonize_trait_names(Trait),
    Model = factor(
      Model,
      levels = c("Population-level", "Between-family", "Within-family")
    ),
    Outcome_disp = Outcome
  )

if (nrow(dat) == 0) {
  stop("No valid rows available for plotting after harmonization.")
}

# ---------- Outcome order ----------
present_outcomes <- preferred_outcome_order[
  preferred_outcome_order %in% unique(dat$Outcome_disp)
]

other_outcomes <- setdiff(
  sort(unique(dat$Outcome_disp)),
  preferred_outcome_order
)

outcome_levels <- rev(c(present_outcomes, other_outcomes))

dat <- dat %>%
  mutate(
    Outcome_disp = factor(Outcome_disp, levels = outcome_levels)
  )

# ---------- Plot function ----------
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
  
  ggplot(df_trait, aes(x = OR, y = Outcome_disp, color = Model, shape = Model)) +
    geom_vline(xintercept = 1, linetype = 2) +
    geom_errorbar(
      aes(
        xmin = CI_low_capped,
        xmax = CI_high_capped
      ),
      orientation = "y",
      position = pd,
      width = 0.25,
      linewidth = 0.4
    )+
    geom_point(position = pd, size = 2.3) +
    scale_x_log10(
      limits = c(lo, hi),
      breaks = c(0.5, 0.67, 1, 1.5, 2, 3),
      labels = scales::label_number(accuracy = 0.01)
    ) +
    scale_shape_manual(
      values = c(
        "Population-level" = 16,
        "Between-family"   = 17,
        "Within-family"    = 4
      )
    ) +
    scale_color_manual(
      values = c(
        "Population-level" = "#F8766D",
        "Between-family"   = "#1F4E79",
        "Within-family"    = "#6BAED6"
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
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.text.y = element_text(size = 11),
      legend.title = element_blank()
    )
}

# ---------- Build individual plots ----------
traits <- sort(unique(dat$Trait))

plots <- setNames(vector("list", length(traits)), traits)

for (tr in traits) {
  df_tr <- dat %>% filter(Trait == tr)
  
  p <- make_trait_plot(df_tr, tr)
  plots[[tr]] <- p
  
  safe_name <- sanitize_filename(tr)
  
  ggsave(
    filename = file.path(out_dir_ind, paste0("forest_EUR_3models_", safe_name, ".png")),
    plot = p,
    width = 6,
    height = 6.5,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(out_dir_ind, paste0("forest_EUR_3models_", safe_name, ".pdf")),
    plot = p,
    width = 6,
    height = 6.5
  )
}

# ---------- Panel pages ----------
chunk_vec <- function(x, n) {
  split(x, ceiling(seq_along(x) / n))
}

page_chunks <- chunk_vec(traits, 8)
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
  
  if (length(page_plots) < 8) {
    page_plots <- c(
      page_plots,
      replicate(
        8 - length(page_plots),
        ggplot() + theme_void(),
        simplify = FALSE
      )
    )
  }
  
  panel_i <- wrap_plots(page_plots, ncol = 2) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  ggsave(
    filename = file.path(out_dir_panel, sprintf("forest_EUR_3models_traits_page_%02d.png", i)),
    plot = panel_i,
    width = 14,
    height = 20,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(out_dir_panel, sprintf("forest_EUR_3models_traits_page_%02d.pdf", i)),
    plot = panel_i,
    width = 14,
    height = 20
  )
  
  panel_list[[i]] <- panel_i
}

# ---------- Combined PDF ----------
pdf_path <- file.path(out_dir_panel, "forest_EUR_3models_ALL.pdf")

pdf(pdf_path, width = 14, height = 20)
for (p in panel_list) {
  print(p)
}
dev.off()

message("Done.")
message("Individual EUR plots saved to: ", out_dir_ind)
message("Panel EUR plots saved to: ", out_dir_panel)
message("Combined EUR PDF saved to: ", pdf_path)
