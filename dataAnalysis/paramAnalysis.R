# ==== Load Libraries ====
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

output_dir <- "./paramPlot/"

# ==== Create Output Folder if Not Exists ====
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

seqLength <- 10000

# ==== Function to fill sites from reference scripts ====
fill_sites_from_reference <- function(df) {
  # Fill wtErrModel from ErrModel
  err_model_sites <- df$sites[df$script == "ErrModel"]
  if(length(err_model_sites) > 0){
    df$sites[df$script == "wtErrModel"] <- err_model_sites[1]
  }
  
  # Fill TruewtErrModel from TrueErrModel
  true_err_sites <- df$sites[df$script == "TrueErrModel"]
  if(length(true_err_sites) > 0){
    df$sites[df$script == "TruewtErrModel"] <- true_err_sites[1]
  }
  
  return(df)
}

analyse_params <- function(model_name, df, param) {
  message("Analysing model: ", model_name)
  
  # ---- Merge parameter combinations ----
  required_param_cols <- c("r_value", "minD_value", "minVar_value")
  if(!all(required_param_cols %in% names(param))) stop("param missing columns")
  
  param_comb <- param[, required_param_cols] %>%
    rename(r_Value = r_value, minD = minD_value, minVar = minVar_value) %>%
    mutate(across(everything(), as.numeric)) %>%
    mutate(param_combo = paste0("D", minD, "V", minVar))
  
  df <- df %>%
    mutate(r_Value = as.numeric(as.character(r_Value))) %>%
    left_join(param_comb, by = "r_Value")
  
  # ---- Scale rates and compute errors ----
  df <- df %>%
    mutate(
      effective_rootRate       = true_rootRate / seqLength,
      effective_tumourRate     = true_tumourRate / seqLength,
      effective_est_rootRate   = est_rootRate / sites,
      effective_est_tumourRate = est_tumourRate / sites,
      error_rootRate           = abs(effective_est_rootRate - effective_rootRate),
      error_tumourRate         = abs(effective_est_tumourRate - effective_tumourRate),
      error_lambda             = abs(est_lambda - true_lambda)
    )
  
  # ---- Summary Stats by Script ----
  stats <- df %>%
    group_by(script) %>%
    summarise(
      mean_error_rootRate   = mean(error_rootRate, na.rm=TRUE),
      rmse_rootRate         = sqrt(mean(error_rootRate^2, na.rm=TRUE)),
      mean_error_tumourRate = mean(error_tumourRate, na.rm=TRUE),
      rmse_tumourRate       = sqrt(mean(error_tumourRate^2, na.rm=TRUE)),
      mean_error_lambda     = mean(error_lambda, na.rm=TRUE),
      rmse_lambda           = sqrt(mean(error_lambda^2, na.rm=TRUE)),
      .groups="drop"
    )
  write.csv(stats, file=file.path(output_dir, paste0(model_name, "_summary_stats_by_script.csv")), row.names=FALSE)
  
  # ---- Scatter Plots by Script ----
  p1 <- ggplot(df, aes(effective_rootRate, effective_est_rootRate, color=script)) +
    geom_point(alpha=0.6) + geom_smooth(aes(group=script), method="lm", se=FALSE) +
    labs(x="True root rate per site", y="Estimated root rate per site", title=paste(model_name, ": Root Rate")) +
    theme_minimal()
  
  p2 <- ggplot(df, aes(effective_tumourRate, effective_est_tumourRate, color=script)) +
    geom_point(alpha=0.6) + geom_smooth(aes(group=script), method="lm", se=FALSE) +
    labs(x="True tumour rate per site", y="Estimated tumour rate per site", title=paste(model_name, ": Tumour Rate")) +
    theme_minimal()
  
  if(model_name != "realistic") {
    p3 <- ggplot(df, aes(true_lambda, est_lambda, color=script)) +
      geom_point(alpha=0.6) + labs(x="True Lambda", y="Estimated Lambda", title=paste(model_name, ": Lambda")) +
      theme_minimal()
    ggsave(file.path(output_dir, paste0(model_name, "_lambda.png")), p3, width=6, height=5, dpi=300)
  }
  
  ggsave(file.path(output_dir, paste0(model_name, "_rootRate.png")), p1, width=6, height=5, dpi=300)
  ggsave(file.path(output_dir, paste0(model_name, "_tumourRate.png")), p2, width=6, height=5, dpi=300)
  
  # ---- Error Distributions ----
  error_plots <- list(
    list(col="error_rootRate", file_suffix="rootRate"),
    list(col="error_tumourRate", file_suffix="tumourRate")
  )
  
  for(e in error_plots){
    p_hist <- ggplot(df, aes_string(x=e$col, fill="script", color="script")) +
      geom_histogram(alpha=0.4, position="identity", bins=40) +
      labs(title=paste(model_name, ":", e$col, "Error Distribution"), x="Absolute Error", y="Count") +
      theme_minimal()
    
    p_dens <- ggplot(df, aes_string(x=e$col, fill="script", color="script")) +
      geom_density(alpha=0.3) +
      labs(title=paste(model_name, ":", e$col, "Error Density"), x="Absolute Error", y="Density") +
      theme_minimal()
    
    ggsave(file.path(output_dir, paste0(model_name, "_error_", e$file_suffix, "_by_script.png")), p_hist, width=5, height=5, dpi=300)
    ggsave(file.path(output_dir, paste0(model_name, "_error_", e$file_suffix, "_density_by_script.png")), p_dens, width=5, height=5, dpi=300)
  }
  
  # ---- Compute mean error per param_combo and script ----
  df_summary <- df %>%
    group_by(script, param_combo, minD, minVar) %>%
    summarise(
      mean_error_rootRate   = mean(error_rootRate, na.rm = TRUE),
      mean_error_tumourRate = mean(error_tumourRate, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(minD, minVar)  # ensure proper order
  
  # Make param_combo a factor in the correct order
  df_summary <- df_summary %>%
    mutate(param_combo = factor(param_combo, levels = unique(param_combo)))
  
  # ---- Plot Root Rate Error ----
  p_root <- ggplot() +
    geom_point(data = df, aes(x = param_combo, y = error_rootRate, color = script), alpha = 0.6) +
    geom_line(data = df_summary, aes(x = param_combo, y = mean_error_rootRate, color = script, group = script), size = 1) +
    labs(title=paste(model_name, ": Root Rate Error by Parameter"),
         x="Parameter Combination", y="Absolute Error") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # ---- Plot Tumour Rate Error ----
  p_tumour <- ggplot() +
    geom_point(data = df, aes(x = param_combo, y = error_tumourRate, color = script), alpha = 0.6) +
    geom_line(data = df_summary, aes(x = param_combo, y = mean_error_tumourRate, color = script, group = script), size = 1) +
    labs(title=paste(model_name, ": Tumour Rate Error by Parameter"),
         x="Parameter Combination", y="Absolute Error") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(output_dir, paste0(model_name, "_error_rootRate_by_param.png")), p_root, width=7, height=5, dpi=300)
  ggsave(file.path(output_dir, paste0(model_name, "_error_tumourRate_by_param.png")), p_tumour, width=7, height=5, dpi=300)
}


# ==== Load Data ====
cat("Loading data files...\n")
simpleData <- read.csv("../simpleModel/4-analysis/summary/parameter_results.csv")
realisticData <- read.csv("../realisticModel/4-analysis/summary/parameter_results.csv")
calibratedData <- read.csv("../calibratedModel/4-analysis/summary/parameter_results.csv")

simpleParam <- read.csv("../simpleModel/4-analysis/summary/var_counts.csv")
realisticParam <- read.csv("../realisticModel/4-analysis/summary/var_counts.csv")
calibratedParam <- read.csv("../calibratedModel/4-analysis/summary/var_counts.csv")

simpleData <- fill_sites_from_reference(simpleData)
realisticData <- fill_sites_from_reference(realisticData)
calibratedData <- fill_sites_from_reference(calibratedData)

# ==== Run Analysis ====
analyse_params("simple", simpleData, simpleParam)
analyse_params("realistic", realisticData, realisticParam)
analyse_params("calibrated", calibratedData, calibratedParam)
