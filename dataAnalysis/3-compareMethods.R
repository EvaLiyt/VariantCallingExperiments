# ===============================
# ---- load libraries ----
# ===============================
library(tidyverse)

# ===============================
# ---- output folder ----
# ===============================
output_dir <- "./compareFilter"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ===============================
# ---- global settings ----
# ===============================
methods <- c("method1", "method2", "method3")

models <- list(
  Simple     = "../simpleModel/4-analysis/summary",
  Calibrated = "../calibratedModel/4-analysis/summary",
  Complex  = "../realisticModel/4-analysis/summary"
)

TOP_N <- 3

# ===============================
# ---- containers ----
# ===============================
best_overall_all_models    <- tibble()
best_per_method_all_models <- tibble()

# ===============================
# ---- helper: load parameter map ----
# ===============================
load_param_map <- function(summary_dir) {
  read.csv(
    file.path(summary_dir, "summary_method1.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) %>%
    select(r_value, minD_value, minVar_value) %>%
    distinct() %>%
    mutate(
      combination = paste0("D", minD_value, "_V", minVar_value)
    ) %>%
    select(r_value, combination)
}

# ===============================
# ---- helper: load F1 ----
# ===============================
load_method_scores <- function(summary_dir, method) {
  read.csv(
    file.path(summary_dir, paste0("summary_", method, ".csv")),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) %>%
    mutate(
      method = method,
      F1_score = Mean_F1Score
    ) %>%
    select(r_value, method, F1_score)
}

# ===============================
# ---- main loop ----
# ===============================
for (model_name in names(models)) {

  message("Processing model: ", model_name)
  summary_dir <- models[[model_name]]

  param_map <- load_param_map(summary_dir)

  raw_scores <- bind_rows(
    lapply(methods, load_method_scores, summary_dir = summary_dir)
  )

  full_data <- raw_scores %>%
    left_join(param_map, by = "r_value")

  # ===============================
  # ---- aggregate mean F1 ----
  # ===============================
  f1_summary <- full_data %>%
    group_by(method, combination) %>%
    summarise(
      mean_F1 = mean(F1_score, na.rm = TRUE),
      sd_F1   = sd(F1_score, na.rm = TRUE),
      n       = sum(!is.na(F1_score)),
      .groups = "drop"
    )

  # ===============================
  # ---- TOP N PER METHOD ----
  # ===============================
  for (m in methods) {

    topN <- f1_summary %>%
      filter(method == m) %>%
      arrange(desc(mean_F1)) %>%
      slice_head(n = TOP_N) %>%
      mutate(
        rank  = row_number(),
        model = model_name
      ) %>%
      select(
        model, method, rank, combination,
        mean_F1, sd_F1, n
      )

    best_per_method_all_models <- bind_rows(
      best_per_method_all_models,
      topN
    )
  }

  # ===============================
  # ---- TOP N OVERALL ----
  # ===============================
  best_overall <- f1_summary %>%
    arrange(desc(mean_F1)) %>%
    slice_head(n = TOP_N) %>%
    mutate(
      rank  = row_number(),
      model = model_name
    ) %>%
    select(
      model, method, rank, combination,
      mean_F1, sd_F1, n
    )

  best_overall_all_models <- bind_rows(
    best_overall_all_models,
    best_overall
  )

  # ===============================
  # ---- plot ----
  # ===============================
  p <- ggplot(
    f1_summary,
    aes(x = combination, y = mean_F1, color = method, group = method)
  ) +
    geom_point(size = 3) +
    geom_line(linewidth = 0.8) +
    labs(
      title = paste0("Mean F1 by Parameter and Method\n(", model_name, " Model)"),
      x = "Parameter Combination",
      y = "Mean F1",
      color = "Method"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  ggsave(
    file.path(output_dir, paste0("F1_parameter_method_", model_name, ".png")),
    p,
    width = 10,
    height = 5,
    dpi = 300
  )
}

# ===============================
# ---- save outputs ----
# ===============================
write_csv(
  best_overall_all_models,
  file.path(output_dir, "TOP3_overall_all_models.csv")
)

write_csv(
  best_per_method_all_models,
  file.path(output_dir, "TOP3_per_method_all_models.csv")
)

# ===============================
# ---- console output ----
# ===============================
cat("\n===== TOP 3 PER METHOD (ALL MODELS) =====\n")
print(best_per_method_all_models)

cat("\nAll models processed successfully.\n")
