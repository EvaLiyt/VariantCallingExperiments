# Install cowplot if not already installed
# if (!require("cowplot", quietly = TRUE)) {
#   cat("Installing cowplot package...\n")
#   install.packages("cowplot", repos = "https://cloud.r-project.org/")
#   cat("cowplot installed successtumoury\n\n")
# }

# Load necessary packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(cowplot)

# Read the results CSV
results_df <- read.csv("./summary/tree_analysis_results_40_clades.csv")

# Define output directory for plots
output_dir <- "./summary/tumourAnalysis"

# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n\n")
}

# Separate True models (no DV) from DV models
results_df <- results_df %>%
  mutate(
    is_true_model = (script == "TrueErrModel" | script == "TruewtErrModel"),
    minD_numeric = as.numeric(minD),
    minVar_numeric = as.numeric(minVar)
  )

# Create x-axis labels
results_df <- results_df %>%
  mutate(
    x_axis_label = ifelse(
      is_true_model,
      script,  # Just use the script name (TrueErrModel or TruewtErrModel)
      paste0("D", minD_numeric, "V", minVar_numeric)  # DV format for others
    )
  )

# Sort the x-axis: TrueErrModel and TruewtErrModel first, then DV combinations sorted by minD and minVar
true_labels <- c("TrueErrModel", "TruewtErrModel")

# Get DV combinations sorted by minD, then minVar
dv_data <- results_df %>%
  filter(!is_true_model) %>%
  distinct(minD_numeric, minVar_numeric, x_axis_label) %>%
  arrange(minD_numeric, minVar_numeric)

dv_labels <- dv_data$x_axis_label

# Create ordered x-axis with all labels
x_axis_labels_ordered <- c(true_labels, dv_labels)
results_df$x_axis_label <- factor(results_df$x_axis_label, levels = x_axis_labels_ordered)

# Get unique models
models <- sort(unique(results_df$model_folder))

# Get unique scripts
scripts <- sort(unique(results_df$script))

# Define colors for scripts
script_colors <- scales::hue_pal()(length(scripts))
names(script_colors) <- scripts

# Create medium versions of script colors for fill (between light scatter and dark)
medium_script_colors <- sapply(script_colors, function(col) {
  rgb_vals <- col2rgb(col)
  # Blend with white at 40% to get a medium shade
  rgb(rgb_vals[1] * 0.7 + 255 * 0.3, 
      rgb_vals[2] * 0.7 + 255 * 0.3, 
      rgb_vals[3] * 0.7 + 255 * 0.3, maxColorValue = 255)
})
names(medium_script_colors) <- scripts

# =======================
# CREATE SUMMARY STATISTICS
# =======================
cat("Generating summary statistics...\n")

summary_stats <- results_df %>%
  group_by(model_folder, x_axis_label, script) %>%
  summarise(
    mean_RF = mean(RF_Dist, na.rm = TRUE),
    sd_RF = sd(RF_Dist, na.rm = TRUE),
    mean_MAE_H = mean(MAE_treeH, na.rm = TRUE),
    sd_MAE_H = sd(MAE_treeH, na.rm = TRUE),
    mean_MAE_L = mean(MAE_treeL, na.rm = TRUE),
    sd_MAE_L = sd(MAE_treeL, na.rm = TRUE),
    mean_true_H = mean(true_treeH, na.rm = TRUE),
    mean_true_L = mean(true_treeL, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(model_folder, x_axis_label, script)

# Save summary statistics
output_file_stats <- file.path(output_dir, "summary_statistics_tumour_clades.csv")
write.csv(summary_stats, output_file_stats, row.names = FALSE)
cat("Saved summary statistics to", output_file_stats, "\n\n")

# =======================
# PLOTTING - 3x3 GRID (MEAN VALUES ONLY)
# =======================

# Metrics to plot
metrics_to_plot <- list(
  list(metric = "RF_Dist", mean_col = "mean_RF", name = "Robinson-Foulds Distance (Tumour Tips)"),
  list(metric = "MAE_treeH", mean_col = "mean_MAE_H", name = "MAE Tree Height (Tumour Tips)"),
  list(metric = "MAE_treeL", mean_col = "mean_MAE_L", name = "MAE Tree Length (Tumour Tips)")
)

# Create a function to generate plots using mean values only for 3x3 grid
create_plot_means_only <- function(data, model_name, metric_info) {
  
  mean_col <- metric_info$mean_col
  metric_name <- metric_info$name
  
  # Filter summary stats for the current model for means
  model_summary <- data %>%
    filter(model_folder == model_name) %>%
    mutate(x_axis_label = factor(x_axis_label, levels = x_axis_labels_ordered))
  
  # Create the plot with only mean points
  p <- ggplot(model_summary, aes(x = x_axis_label, y = !!sym(mean_col))) +
    # Add mean points (border and fill same medium color)
    geom_point(aes(color = script, fill = script),
               size = 6, shape = 21, stroke = 1, alpha = 0.9) +
    
    # Add vertical line to separate True models from DV models
    geom_vline(xintercept = 2.5, color = "gray50", linetype = "dotted", linewidth = 0.8, alpha = 0.7) +
    
    scale_color_manual(values = medium_script_colors) +
    scale_fill_manual(values = medium_script_colors) +
    
    labs(
      title = metric_name,
      x = "",
      y = metric_name,
      color = "Script",
      fill = "Script"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
      axis.title.y = element_text(size = 9),
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
      legend.position = "none",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.3, "cm"),
      panel.grid.major.y = element_line(color = "lightgray", linewidth = 0.3),
      panel.grid.minor.y = element_line(color = "lightgray", linewidth = 0.15, linetype = "dotted"),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    )
  
  return(p)
}

# Create all plots first
cat("Creating plots for 3x3 grid...\n")

all_plots <- list()
plot_index <- 1

for (metric_idx in seq_along(metrics_to_plot)) {
  for (model_idx in seq_along(models)) {
    
    model_name <- models[model_idx]
    metric_info <- metrics_to_plot[[metric_idx]]
    
    p <- create_plot_means_only(summary_stats, model_name, metric_info)
    all_plots[[plot_index]] <- p
    
    plot_index <- plot_index + 1
  }
}

# Create legend
legend_data <- summary_stats %>%
  distinct(script) %>%
  mutate(fill = script)

dummy_plot <- ggplot(legend_data, aes(x = script, fill = script)) +
  geom_bar() +
  scale_fill_manual(values = medium_script_colors) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  guides(fill = guide_legend(nrow = 1))

legend <- get_legend(dummy_plot)

# Arrange in 3x3 grid (3 rows for metrics, 3 columns for models)
combined_plot <- plot_grid(
  plotlist = all_plots,
  ncol = 3,  # 3 columns for 3 models
  nrow = 3,  # 3 rows for 3 metrics
  align = "hv",
  axis = "both"
)

# Add legend at the bottom
final_plot <- plot_grid(
  combined_plot,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.08)
)

# Save the combined 3x3 plot
output_filename <- file.path(output_dir, "plot_models_tumour_3x3.png")
ggsave(output_filename, final_plot, width = 16, height = 14, dpi = 300)

cat("Saved 3x3 combined plot to", output_filename, "\n\n")

# =======================
# PLOTTING - SCATTER + MEAN OVERLAY FOR COMPARISON PLOTS
# =======================

cat("Creating combined comparison plots...\n")

# Create a combined comparison plot - all models together for RF distance
p_combined_rf <- ggplot(results_df, aes(x = x_axis_label, y = RF_Dist)) +
  geom_jitter(aes(color = script), 
              width = 0.15, 
              alpha = 0.4, 
              size = 2,
              stroke = 0) +
  geom_point(data = summary_stats, aes(x = x_axis_label, y = mean_RF, color = script, fill = script),
             size = 4, shape = 21, stroke = 1, alpha = 0.9) +
  geom_vline(xintercept = 2.5, color = "gray50", linetype = "dotted", linewidth = 0.8, alpha = 0.7) +
  facet_wrap(~factor(model_folder, levels = models), ncol = 3) +
  scale_color_manual(values = medium_script_colors) +
  scale_fill_manual(values = medium_script_colors) +
  labs(
    title = "Robinson-Foulds Distance (Tumour Tips) Across All Models",
    x = "Script / DV Combination",
    y = "Robinson-Foulds Distance",
    color = "Script",
    fill = "Script"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    legend.position = "bottom",
    panel.grid.major.y = element_line(color = "lightgray", linewidth = 0.3),
    strip.text = element_text(face = "bold", size = 10),
    plot.margin = margin(10, 10, 10, 10)
  )

# Save the combined comparison plot for RF distance
output_file_rf <- file.path(output_dir, "plot_combined_comparison_RF_tumour.png")
ggsave(output_file_rf, p_combined_rf, width = 15, height = 6, dpi = 300)
cat("Saved combined comparison RF plot to", output_file_rf, "\n")

# Create a combined comparison plot for MAE Tree Height
p_combined_mae_h <- ggplot(results_df, aes(x = x_axis_label, y = MAE_treeH)) +
  geom_jitter(aes(color = script), 
              width = 0.15, 
              alpha = 0.4, 
              size = 2,
              stroke = 0) +
  geom_point(data = summary_stats, aes(x = x_axis_label, y = mean_MAE_H, color = script, fill = script),
             size = 4, shape = 21, stroke = 1, alpha = 0.9) +
  geom_vline(xintercept = 2.5, color = "gray50", linetype = "dotted", linewidth = 0.8, alpha = 0.7) +
  facet_wrap(~factor(model_folder, levels = models), ncol = 3) +
  scale_color_manual(values = medium_script_colors) +
  scale_fill_manual(values = medium_script_colors) +
  labs(
    title = "MAE Tree Height (Tumour Tips) Across All Models",
    x = "Script / DV Combination",
    y = "MAE Tree Height",
    color = "Script",
    fill = "Script"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    legend.position = "bottom",
    panel.grid.major.y = element_line(color = "lightgray", linewidth = 0.3),
    strip.text = element_text(face = "bold", size = 10),
    plot.margin = margin(10, 10, 10, 10)
  )

# Save the MAE tree height comparison plot
output_file_mae_h <- file.path(output_dir, "plot_combined_comparison_MAE_treeH_tumour.png")
ggsave(output_file_mae_h, p_combined_mae_h, width = 15, height = 6, dpi = 300)
cat("Saved MAE tree height comparison plot to", output_file_mae_h, "\n")

# Create a combined comparison plot for MAE Tree Length
p_combined_mae_l <- ggplot(results_df, aes(x = x_axis_label, y = MAE_treeL)) +
  geom_jitter(aes(color = script), 
              width = 0.15, 
              alpha = 0.4, 
              size = 2,
              stroke = 0) +
  geom_point(data = summary_stats, aes(x = x_axis_label, y = mean_MAE_L, color = script, fill = script),
             size = 4, shape = 21, stroke = 1, alpha = 0.9) +
  geom_vline(xintercept = 2.5, color = "gray50", linetype = "dotted", linewidth = 0.8, alpha = 0.7) +
  facet_wrap(~factor(model_folder, levels = models), ncol = 3) +
  scale_color_manual(values = medium_script_colors) +
  scale_fill_manual(values = medium_script_colors) +
  labs(
    title = "MAE Tree Length (Tumour Tips) Across All Models",
    x = "Script / DV Combination",
    y = "MAE Tree Length",
    color = "Script",
    fill = "Script"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    legend.position = "bottom",
    panel.grid.major.y = element_line(color = "lightgray", linewidth = 0.3),
    strip.text = element_text(face = "bold", size = 10),
    plot.margin = margin(10, 10, 10, 10)
  )

# Save the MAE tree length comparison plot
output_file_mae_l <- file.path(output_dir, "plot_combined_comparison_MAE_treeL_tumour.png")
ggsave(output_file_mae_l, p_combined_mae_l, width = 15, height = 6, dpi = 300)
cat("Saved MAE tree length comparison plot to", output_file_mae_l, "\n")

cat("\n===========================================\n")
cat("PLOTTING COMPLETED SUCCESStumourY!\n")
cat("===========================================\n\n")
cat("Generated files saved to:", output_dir, "\n\n")
cat("Main plot (3x3 grid - mean values only):\n")
cat("- plot_models_tumour_3x3.png (3 metrics x 3 models in one grid)\n\n")
cat("Combined comparison plots (scatter + mean overlay with facets):\n")
cat("- plot_combined_comparison_RF_tumour.png\n")
cat("- plot_combined_comparison_MAE_treeH_tumour.png\n")
cat("- plot_combined_comparison_MAE_treeL_tumour.png\n\n")
cat("Summary statistics:\n")
cat("- summary_statistics_tumour_clades.csv\n\n")
cat("Features:\n")
cat("- 3x3 grid: shows only mean values (no scatter)\n")
cat("- Mean dots: medium shade (between light scatter and original color)\n")
cat("- Border and fill: same medium color\n")
cat("- Comparison plots: scatter + mean overlay\n")
cat("- Single legend at bottom\n")