# Install cowplot if not already installed
# if (!require("cowplot", quietly = TRUE)) {
#   cat("Installing cowplot package...\n")
#   install.packages("cowplot", repos = "https://cloud.r-project.org/")
#   cat("cowplot installed successfully\n\n")
# }

# Load necessary packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(cowplot)

# Create output directory if it doesn't exist
dir.create("./summary/param", showWarnings = FALSE, recursive = TRUE)

# Read the parameter analysis results
param_results <- read.csv("./summary/parameter_analysis_results.csv")

# Remove rows with NA DV values for non-True models
param_results <- param_results %>%
  mutate(
    is_true_model = (script == "TrueErrModel" | script == "TruewtErrModel"),
    minD_numeric = as.numeric(minD),
    minVar_numeric = as.numeric(minVar)
  )

# Create x-axis labels
param_results <- param_results %>%
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
dv_data <- param_results %>%
  filter(!is_true_model) %>%
  distinct(minD_numeric, minVar_numeric, x_axis_label) %>%
  arrange(minD_numeric, minVar_numeric)

dv_labels <- dv_data$x_axis_label

# Create ordered x-axis with all labels
x_axis_labels_ordered <- c(true_labels, dv_labels)
param_results$x_axis_label <- factor(param_results$x_axis_label, levels = x_axis_labels_ordered)

# Get unique models
models <- sort(unique(param_results$model_folder))

# Get unique scripts
scripts <- sort(unique(param_results$script))

# Define colors for scripts
script_colors <- scales::hue_pal()(length(scripts))
names(script_colors) <- scripts

# Create medium versions of script colors (between light scatter and original)
medium_script_colors <- sapply(script_colors, function(col) {
  rgb_vals <- col2rgb(col)
  # Blend with white at 40% to get a medium shade
  rgb(rgb_vals[1] * 0.7 + 255 * 0.3, 
      rgb_vals[2] * 0.7 + 255 * 0.3, 
      rgb_vals[3] * 0.7 + 255 * 0.3, maxColorValue = 255)
})
names(medium_script_colors) <- scripts

# Parameters to plot (estimated vs true)
params_to_plot <- list(
  list(est = "est_lambda", true = "true_lambda", error = "error_lambda", 
       name = "Lambda", error_name = "Lambda Error"),
  list(est = "est_rootRate", true = "true_rootRate", error = "error_rootRate", 
       name = "Root Rate", error_name = "Root Rate Error"),
  list(est = "est_tumourRate", true = "true_tumourRate", error = "error_tumourRate", 
       name = "Tumour Rate", error_name = "Tumour Rate Error"),
  list(est = "est_mut_rate_per_site", true = "true_mut_rate_per_site", error = "error_mut_rate_per_site", 
       name = "Mutation Rate per Site", error_name = "Mutation Rate per Site Error")
)

# Create a function to generate plots for estimated vs true values
create_comparison_plot <- function(data, model_name, param_info) {
  
  est_col <- param_info$est
  true_col <- param_info$true
  param_name <- param_info$name
  
  # Filter data for the current model
  model_data <- data %>%
    filter(model_folder == model_name) %>%
    filter(!is.na(!!sym(est_col)) & !is.na(!!sym(true_col)))
  
  # Create the plot
  p <- ggplot(model_data, aes(x = x_axis_label)) +
    # Plot true values as horizontal lines per DV combination
    geom_jitter(aes(y = !!sym(true_col), color = "True"), 
                width = 0.1, 
                alpha = 0.3, 
                size = 2,
                stroke = 0,
                shape = 4) +
    # Plot estimated values as points
    geom_jitter(aes(y = !!sym(est_col), color = script), 
                width = 0.2, 
                alpha = 0.4, 
                size = 2.5,
                stroke = 0) +
    # Add mean points with medium script color
    stat_summary(aes(y = !!sym(est_col), color = script, fill = script), 
                 fun = mean, 
                 geom = "point", 
                 size = 6, 
                 shape = 21,
                 stroke = 1,
                 alpha = 0.9) +
    
    # Add vertical line to separate True models from DV models
    geom_vline(xintercept = 2.5, color = "gray50", linetype = "dotted", linewidth = 0.8, alpha = 0.7) +
    
    scale_color_manual(values = c(script_colors, "True" = "red")) +
    scale_fill_manual(values = medium_script_colors, guide = "none") +
    
    labs(
      title = paste(param_name, "- Estimated vs True"),
      x = "",
      y = param_name,
      color = "Source"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
      axis.title.y = element_text(size = 9),
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
      legend.position = "bottom",
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

# Create a function to generate error plots
create_error_plot <- function(data, model_name, param_info) {
  
  error_col <- param_info$error
  param_name <- param_info$error_name
  
  # Filter data for the current model
  model_data <- data %>%
    filter(model_folder == model_name) %>%
    filter(!is.na(!!sym(error_col)))
  
  # Create the plot
  p <- ggplot(model_data, aes(x = x_axis_label, y = !!sym(error_col))) +
    # Add scatter points with different colors for scripts
    geom_jitter(aes(color = script), 
                width = 0.2, 
                alpha = 0.4, 
                size = 2.5,
                stroke = 0) +
    # Add mean points with medium script color
    stat_summary(aes(color = script, fill = script),
                 fun = mean, 
                 geom = "point", 
                 size = 6, 
                 shape = 21,
                 stroke = 1,
                 alpha = 0.9) +
    # Add zero line reference
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 0.9, alpha = 0.5) +
    
    # Add vertical line to separate True models from DV models
    geom_vline(xintercept = 2.5, color = "gray50", linetype = "dotted", linewidth = 0.8, alpha = 0.7) +
    
    scale_color_manual(values = script_colors) +
    scale_fill_manual(values = medium_script_colors, guide = "none") +
    
    labs(
      title = param_name,
      x = "",
      y = param_name,
      color = "Script"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
      axis.title.y = element_text(size = 9),
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
      legend.position = "bottom",
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

# Create plots for each model
for (model_name in models) {
  
  cat("Creating parameter plots for", model_name, "...\n")
  
  # Create comparison plots (estimated vs true)
  comparison_plots <- list()
  for (i in seq_along(params_to_plot)) {
    p <- create_comparison_plot(param_results, model_name, params_to_plot[[i]])
    comparison_plots[[i]] <- p
  }
  
  # Arrange the three comparison plots
  comparison_combined <- plot_grid(
    plotlist = comparison_plots,
    ncol = 3,
    align = "h",
    axis = "b"
  )
  
  # Add title for the comparison plots
  comparison_final <- plot_grid(
    comparison_combined,
    top = paste("Model:", model_name, "- Estimated vs True Values"),
    rel_heights = c(1, 0.05)
  )
  
  # Save the comparison plot
  output_filename <- paste0("./summary/param/param_comparison_", gsub(" ", "_", model_name), ".png")
  ggsave(output_filename, comparison_final, width = 16, height = 5, dpi = 300)
  cat("  Saved comparison plot to", output_filename, "\n")
  
  # Create error plots
  error_plots <- list()
  for (i in seq_along(params_to_plot)) {
    p <- create_error_plot(param_results, model_name, params_to_plot[[i]])
    error_plots[[i]] <- p
  }
  
  # Arrange the three error plots
  error_combined <- plot_grid(
    plotlist = error_plots,
    ncol = 3,
    align = "h",
    axis = "b"
  )
  
  # Add title for the error plots
  error_final <- plot_grid(
    error_combined,
    top = paste("Model:", model_name, "- Estimation Errors"),
    rel_heights = c(1, 0.05)
  )
  
  # Save the error plot
  error_filename <- paste0("./summary/param/param_errors_", gsub(" ", "_", model_name), ".png")
  ggsave(error_filename, error_final, width = 16, height = 5, dpi = 300)
  cat("  Saved error plot to", error_filename, "\n")
}

# Create combined comparison plots - all models together
for (i in seq_along(params_to_plot)) {
  
  param_info <- params_to_plot[[i]]
  est_col <- param_info$est
  true_col <- param_info$true
  param_name <- param_info$name
  
  # Prepare data
  plot_data <- param_results %>%
    filter(!is.na(!!sym(est_col)) & !is.na(!!sym(true_col)))
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = x_axis_label)) +
    geom_jitter(aes(y = !!sym(est_col), color = script), 
                width = 0.15, 
                alpha = 0.4, 
                size = 2,
                stroke = 0) +
    stat_summary(aes(y = !!sym(est_col), color = script, fill = script),
                 fun = mean, 
                 geom = "point", 
                 size = 5, 
                 shape = 21,
                 stroke = 1,
                 alpha = 0.9) +
    geom_vline(xintercept = 2.5, color = "gray50", linetype = "dotted", linewidth = 0.8, alpha = 0.7) +
    facet_wrap(~factor(model_folder, levels = models), ncol = 3) +
    scale_color_manual(values = script_colors) +
    scale_fill_manual(values = medium_script_colors, guide = "none") +
    labs(
      title = paste(param_name, "- Estimated Values Across All Models"),
      x = "Script / DV Combination",
      y = param_name,
      color = "Script"
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
  
  # Save the combined comparison plot
  combined_filename <- paste0("./summary/param/param_combined_", tolower(gsub(" ", "_", param_name)), ".png")
  ggsave(combined_filename, p, width = 15, height = 6, dpi = 300)
  cat("Saved combined", param_name, "plot to", combined_filename, "\n")
}

cat("\n===========================================\n")
cat("PARAMETER PLOTTING COMPLETED\n")
cat("===========================================\n")
cat("\nGenerated plots:\n")
cat("- param_comparison_simpleModel.png\n")
cat("- param_comparison_complexModel.png\n")
cat("- param_comparison_CalibratedModel.png\n")
cat("- param_errors_simpleModel.png\n")
cat("- param_errors_complexModel.png\n")
cat("- param_errors_CalibratedModel.png\n")
cat("- param_combined_lambda.png\n")
cat("- param_combined_root_rate.png\n")
cat("- param_combined_tumour_rate.png\n")
cat("- param_combined_mutation_rate_per_site.png\n")
cat("\nAll plots saved to ./summary/param/\n")
cat("\nFeatures:\n")
cat("- Mean points use medium shade (between light scatter and original color)\n")
cat("- Border and fill same color for mean points\n")