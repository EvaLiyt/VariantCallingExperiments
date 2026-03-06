# ---- install packages if needed ----
# install.packages(c("tidyverse", "viridis", "reshape2", "ggpmisc", "patchwork", "pROC", "RColorBrewer", "ggrepel"))

# ---- load libraries ----
library(tidyverse)
library(dplyr)
library(viridis)
library(reshape2)
library(ggpmisc)
library(patchwork)
library(pROC)
library(RColorBrewer)
library(ggrepel)

# ---- output folders ----
output_dir <- "./performancePlot"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

candidatePlot_dir <- "./candidatePlot"
if (!dir.exists(candidatePlot_dir)) dir.create(candidatePlot_dir, recursive = TRUE)

# ---- swap tumour/normal columns (kept; not used below unless you call it) ----
swap_tumour_normal_csv <- function(filepath) {
  df <- read.csv(filepath, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(df) <- gsub("tumour", "TEMP_TUMOUR", colnames(df))
  colnames(df) <- gsub("normal", "tumour", colnames(df))
  colnames(df) <- gsub("TEMP_TUMOUR", "normal", colnames(df))
  df[] <- lapply(df, function(col) {
    if (is.character(col)) {
      col <- gsub("tumour", "TEMP_TUMOUR", col)
      col <- gsub("normal", "tumour", col)
      col <- gsub("TEMP_TUMOUR", "normal", col)
    }
    return(col)
  })
  write.csv(df, file = filepath, row.names = FALSE, quote = TRUE)
}
# swap_tumour_normal_csv("../simpleModel/4-analysis/summary/var_counts_method2.csv")

# ---- helper: normalize cell names for tumour/normal logic ----
# remove leading r<number>_ so r0_clade_11 -> clade_11 ; r12_cellA -> cellA
clean_cell_name <- function(x) {
  x <- as.character(x)
  sub("^r\\d+_", "", x)
}

# tumour if contains "clade_" after removing rX_
is_tumour <- function(x) {
  x <- clean_cell_name(x)
  grepl("clade_", x)
}

# ---- analyse function ----
analyse_model <- function(model_name, df, vars, hammingDist) {

  # --- ensure numeric + combination ---
  df <- df %>%
    mutate(
      minD_value   = as.numeric(minD_value),
      minVar_value = as.numeric(minVar_value),
      combination  = paste0("D", minD_value, "_V", minVar_value),
      F1_score     = Mean_F1Score   # unify: always use Python F1
    ) %>%
    arrange(minD_value, minVar_value)

  # --- summary per parameter ---
  model_summary <- df %>%
    group_by(minD_value, minVar_value, combination) %>%
    summarise(
      avg_tpr       = mean(Mean_TPR, na.rm = TRUE),
      avg_fpr       = mean(Mean_FPR, na.rm = TRUE),
      avg_precision = mean(Mean_Precision, na.rm = TRUE),
      F1_score      = mean(F1_score, na.rm = TRUE),
      count         = n(),
      .groups       = "drop"
    ) %>%
    arrange(minD_value, minVar_value)

  model_summary$combination <- factor(model_summary$combination, levels = unique(model_summary$combination))

  # ---- Plot1: Precision vs TPR ----
  precision_tpr <- ggplot(model_summary,
                          aes(x = avg_tpr, y = avg_precision, color = combination)) +
    geom_point(size = 2, alpha = 0.9) +
    geom_text(aes(label = combination), hjust = 0.5, vjust = -0.7, size = 3) +
    scale_color_viridis_d(option = "C") +
    coord_cartesian(xlim = c(0.2, 0.45), ylim = c(0.55, 1)) +
    labs(
      title = "Precision and True Positive Rate (TPR)",
      subtitle = paste0("For ", model_name, " model. Each point = (minD, minVar) combination"),
      x = "True Positive Rate (TPR)",
      y = "Precision"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      legend.position = "none"
    )

  ggsave(file.path(output_dir, paste0("precision_tpr_", model_name, ".png")),
         plot = precision_tpr, width = 7, height = 6, dpi = 300)

  # ---- Plot 2: F1 Score (raw points + linear fit) ----
  df$combination <- factor(df$combination, levels = levels(model_summary$combination))

  lm_f1 <- lm(F1_score ~ minD_value + minVar_value, data = df)
  df$predicted_F1 <- predict(lm_f1, newdata = df)

  reg_summary <- paste0(
    "y = ", round(coef(lm_f1)[1], 3), " + ",
    round(coef(lm_f1)[2], 3), "·minD + ",
    round(coef(lm_f1)[3], 3), "·minVar\n",
    "R² = ", round(summary(lm_f1)$r.squared, 3)
  )

  F1_scatter_plot <- ggplot(df, aes(x = combination, y = F1_score)) +
    geom_point(color = "#FC8D62", size = 2, alpha = 0.8,
               position = position_jitter(width = 0.3)) +
    annotate("text", x = Inf, y = Inf, label = reg_summary,
             hjust = 1.1, vjust = 1.1, size = 4, color = "black") +
    theme_minimal() +
    labs(
      title = paste0("F1 Score for ", model_name, " Model"),
      subtitle = "F1 from Python (Mean_F1Score), per replicate",
      x = "Parameter Combination (minD_minVar)",
      y = "F1 Score"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(output_dir, paste0("F1_score_linear_fit_", model_name, ".png")),
         plot = F1_scatter_plot, width = 6, height = 5, dpi = 300)

  # ---- hammingDist plot ----
  ham <- hammingDist %>%
    mutate(
      cell1_clean = clean_cell_name(cell1),
      cell2_clean = clean_cell_name(cell2),
      cell1_type  = ifelse(is_tumour(cell1_clean), "C", "N"),
      cell2_type  = ifelse(is_tumour(cell2_clean), "C", "N"),
      pair_type   = case_when(
        cell1_type == "N" & cell2_type == "N" ~ "NN",
        cell1_type == "C" & cell2_type == "C" ~ "CC",
        TRUE ~ "NC"
      )
    )

  hist_overall <- ggplot(ham, aes(x = hammingDist, fill = pair_type)) +
    geom_histogram(alpha = 0.6, position = "identity", bins = 50) +
    theme_minimal() +
    labs(
      title = paste0("Distribution of pairwise hamming distance - ", model_name),
      x = "pairwise hamming distance",
      y = "Count",
      fill = "Pair type"
    )

  ggsave(
    file.path(output_dir, paste0("hamming_distance_", model_name, ".png")),
    plot = hist_overall, width = 7, height = 5, dpi = 300
  )

  df$model <- model_name
  model_summary$model <- model_name

  return(list(summary = model_summary, raw = df))
}

# ---- Load Data ----
simpleData <- read.csv("../simpleModel/4-analysis/summary/summary_method2.csv")
simpleVar <- read.csv("../simpleModel/4-analysis/summary/var_counts_method2.csv")
simplehammingDist <- read.table("../simpleModel/4-analysis/summary/hammingDist_method2.tsv", sep = "\t", header = TRUE)

realisticData <- read.csv("../realisticModel/4-analysis/summary/summary_method2.csv")
realisticVar <- read.csv("../realisticModel/4-analysis/summary/var_counts_method2.csv")
realistichammingDist <- read.table("../realisticModel/4-analysis/summary/hammingDist_method2.tsv", sep = "\t", header = TRUE)

calibratedData <- read.csv("../calibratedModel/4-analysis/summary/summary_method2.csv")
calibratedVar <- read.csv("../calibratedModel/4-analysis/summary/var_counts_method2.csv")
calibratedhammingDist <- read.table("../calibratedModel/4-analysis/summary/hammingDist_method2.tsv", sep = "\t", header = TRUE)

# ---- analyse each model ----
simple_results     <- analyse_model("Simple", simpleData, simpleVar, simplehammingDist)
realistic_results  <- analyse_model("Complex", realisticData, realisticVar, realistichammingDist)
calibrated_results <- analyse_model("Calibrated", calibratedData, calibratedVar, calibratedhammingDist)

simple_summary     <- simple_results$summary
realistic_summary  <- realistic_results$summary
calibrated_summary <- calibrated_results$summary

simple_raw         <- simple_results$raw
realistic_raw      <- realistic_results$raw
calibrated_raw     <- calibrated_results$raw

combined_raw <- bind_rows(simple_raw, realistic_raw, calibrated_raw) %>%
  mutate(model = factor(model, levels = c("Simple", "Calibrated", "Complex")))

combined_summary <- bind_rows(simple_summary, calibrated_summary, realistic_summary) %>%
  mutate(model = factor(model, levels = c("Simple", "Calibrated", "Complex")))

# ---- IMPORTANT: align combination levels across models ----
all_combos <- sort(unique(c(as.character(combined_raw$combination),
                            as.character(combined_summary$combination))))
combined_raw$combination <- factor(as.character(combined_raw$combination), levels = all_combos)
combined_summary$combination <- factor(as.character(combined_summary$combination), levels = all_combos)

# ---- compute regression labels safely (no cur_data/rowwise/list-column pitfalls) ----
regression_labels <- combined_raw %>%
  group_by(model) %>%
  mutate(x = as.numeric(factor(combination, levels = all_combos))) %>%
  summarise(
    intercept = {
      fit <- lm(F1_score ~ x)
      unname(coef(fit)[1])
    },
    slope = {
      fit <- lm(F1_score ~ x)
      unname(coef(fit)[2])
    },
    r2 = {
      fit <- lm(F1_score ~ x)
      summary(fit)$r.squared
    },
    .groups = "drop"
  ) %>%
  mutate(
    formula_label = paste0(
      model, ": y = ", round(intercept, 3),
      " + ", round(slope, 3),
      "·x, R² = ", round(r2, 3)
    )
  ) %>%
  arrange(model) %>%                    # keep stable order
  mutate(y_pos = 1 - 0.08 * (row_number() - 1))  # stack vertically

# ---- combined F1 plot (NO LINES; regression labels stacked vertically) ----
combined_f1_plot <- ggplot() +
  geom_point(
    data = combined_raw,
    aes(x = combination, y = F1_score, color = model),
    alpha = 0.35, size = 1.6,
    position = position_jitter(width = 0.18, height = 0)
  ) +
  geom_point(
    data = combined_summary,
    aes(x = combination, y = F1_score, color = model),
    size = 3
  ) +
  # geom_text(
  #   data = regression_labels,
  #   aes(x = Inf, y = y_pos, label = formula_label, color = model),
  #   hjust = 1.05, vjust = 1, size = 4,
  #   inherit.aes = FALSE,
  #   show.legend = FALSE
  # ) +
  labs(
    title = "F1 Score Comparison Across Models",
    x = "minD_minVar Combination",
    y = "F1 Score"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank()
  )

ggsave(
  file.path(candidatePlot_dir, "f1_overall.png"),
  plot = combined_f1_plot, width = 9, height = 5, dpi = 300
)

# ---- Summary barplot across datasets ----
simple_vars_means <- simpleVar %>%
  summarise(
    avg_called_sites_normal = mean(avg_total_variants_normal, na.rm = TRUE),
    avg_called_var_normal   = mean(avg_non_zero_normal, na.rm = TRUE),
    avg_called_sites_tumour = mean(avg_total_variants_tumour, na.rm = TRUE),
    avg_called_var_tumour   = mean(avg_non_zero_tumour, na.rm = TRUE)
  ) %>% mutate(dataset = "Simple")

realistic_vars_means <- realisticVar %>%
  summarise(
    avg_called_sites_normal = mean(avg_total_variants_normal, na.rm = TRUE),
    avg_called_var_normal   = mean(avg_non_zero_normal, na.rm = TRUE),
    avg_called_sites_tumour = mean(avg_total_variants_tumour, na.rm = TRUE),
    avg_called_var_tumour   = mean(avg_non_zero_tumour, na.rm = TRUE)
  ) %>% mutate(dataset = "Complex")

calibrated_vars_means <- calibratedVar %>%
  summarise(
    avg_called_sites_normal = mean(avg_total_variants_normal, na.rm = TRUE),
    avg_called_var_normal   = mean(avg_non_zero_normal, na.rm = TRUE),
    avg_called_sites_tumour = mean(avg_total_variants_tumour, na.rm = TRUE),
    avg_called_var_tumour   = mean(avg_non_zero_tumour, na.rm = TRUE)
  ) %>% mutate(dataset = "Calibrated")

combined_means <- bind_rows(simple_vars_means, realistic_vars_means, calibrated_vars_means)

combined_long <- combined_means %>%
  pivot_longer(cols = starts_with("avg_"), names_to = "variable", values_to = "value")

summary_var_counts <- ggplot(combined_long, aes(x = variable, y = value, fill = dataset)) +
  geom_col(position = "dodge") +
  labs(
    title = "Comparison of Mean Variant Metrics Across Datasets",
    x = "Metric", y = "Mean Value", fill = "Dataset"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(candidatePlot_dir, "varCounts.png"),
  plot = summary_var_counts, width = 7, height = 5, dpi = 300
)
