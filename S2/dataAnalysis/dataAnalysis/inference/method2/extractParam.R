#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(parallel)

# Install packages if needed
if (!require("data.table", quietly = TRUE)) {
  install.packages("data.table")
  library(data.table)
}

# Define the model folder names
model_folders <- c("simpleModel", "complexModel", "calibratedModel")
base_path <- "../../../"

output_dir <- "./summary"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Model types within each folder
model_types <- c("ErrModel", "wtErrModel", "TrueErrModel", "TruewtErrModel")

# Function to extract true parameter values from simulation log
extract_true_values <- function(filepath) {
  lines <- readLines(filepath, warn = FALSE)
  lines <- lines[grepl("\\S", lines) & !grepl("^#", lines)]
  if (length(lines) < 2) {
    warning(paste("File format error:", filepath))
    return(list(tumourAge=NA, lambda=NA, rootRate=NA, tumourRate=NA, rootAge=NA))
  }
  delimiter <- if (grepl("\t", lines[1])) "\t" else "[[:space:]]+"
  params <- unlist(strsplit(trimws(lines[1]), delimiter))
  values <- as.numeric(unlist(strsplit(trimws(lines[2]), delimiter)))
  if (length(params) != length(values)) {
    warning(paste("Mismatch in", filepath))
    return(list(tumourAge=NA, lambda=NA, rootRate=NA, tumourRate=NA, rootAge=NA))
  }
  get_param_value <- function(prefix) {
    idx <- grep(prefix, params, fixed=TRUE)
    if (length(idx) > 0) values[idx[1]] else NA
  }
  lambda     <- get_param_value("lambda")
  tumourAge  <- get_param_value("tumourAge")
  rootAge    <- get_param_value("rootAge")
  rootRate   <- get_param_value("rootRate")
  tumourRate <- get_param_value("tumourRate")
  if (is.na(tumourRate)) tumourRate <- get_param_value("mu")
  list(tumourAge=tumourAge, lambda=lambda, rootRate=rootRate, tumourRate=tumourRate, rootAge=rootAge)
}

# Function to extract estimated parameter values from result stats file
extract_result_values <- function(filepath) {
  result_data <- read.table(filepath, header = TRUE, stringsAsFactors = FALSE)
  params <- colnames(result_data)
  means <- as.numeric(result_data[1, ])
  get_param_value <- function(name) if (name %in% params) means[params == name] else NA
  lambda     <- get_param_value("lambda")
  rootRate   <- get_param_value("rootRate")
  tree_height <- get_param_value("tree.height")
  tumourRate <- get_param_value("tumourRate")
  if (is.na(tumourRate)) tumourRate <- get_param_value("mu")
  list(lambda=lambda, rootRate=rootRate, tumourRate=tumourRate, tree_height=tree_height)
}

# Function to parse parameters from pileup folder
parse_params <- function(file_name) {
  minD <- NA
  minVar <- NA
  
  if (file.exists(file_name)) {
    lines <- readLines(file_name)
    for (line in lines) {
      line <- trimws(line)
      if (startsWith(line, "min_depth:")) {
        minD <- trimws(sub("min_depth:", "", line))
      } else if (startsWith(line, "min_var:")) {
        minVar <- trimws(sub("min_var:", "", line))
      }
    }
  }
  
  return(list(minD = minD, minVar = minVar))
}

# Function to extract variable sites from err file
extract_variable_sites <- function(err_filepath) {
  if (!file.exists(err_filepath)) {
    return(NA)
  }
  
  lines <- readLines(err_filepath, warn = FALSE)
  
  # Find the first line starting with "INFO: Extract"
  info_lines <- grep("^INFO: Extract", lines, value = TRUE)
  
  if (length(info_lines) == 0) {
    return(NA)
  }
  
  # Extract the first matching line
  info_line <- info_lines[1]
  
  # Parse the line to extract variable sites
  # Pattern: "INFO: Extract XXX variable sites from YYYY sites in alignment RA"
  match <- regmatches(info_line, regexpr("Extract ([0-9]+) variable sites", info_line))
  
  if (length(match) > 0) {
    # Extract the number
    var_sites <- as.integer(sub("Extract ([0-9]+) variable sites", "\\1", match[1]))
    return(var_sites)
  }
  
  return(NA)
}

# Process a single r value and model
process_r_model <- function(r_value, model, model_folder) {
  trueLogFolder <- file.path(base_path, model_folder, "1-simulation/simulationData")
  resultFolder  <- file.path(base_path, model_folder, "3-xml-method2/parsed_results")
  paramsFolder  <- file.path(base_path, model_folder, "1-simulation/pileup")
  xmlScriptFolder <- file.path(base_path, model_folder, "3-xml-method2", paste0(model, "Scripts"))
  
  true_log_file <- file.path(trueLogFolder, paste0("sequenceGeneration_r", r_value, ".log"))
  result_file   <- file.path(resultFolder, paste0(model, "_r", r_value, "_stats.log"))
  param_file    <- file.path(paramsFolder, paste0("r", r_value), paste0("parameters_r", r_value, ".txt"))
  err_file      <- file.path(xmlScriptFolder, paste0(model, "_r", r_value, ".err"))
  
  if (file.exists(true_log_file) && file.exists(result_file)) {
    true_values   <- extract_true_values(true_log_file)
    result_values <- extract_result_values(result_file)
    params        <- parse_params(param_file)
    
    # Extract variable sites and calculate mutation rates per site
    variable_sites <- extract_variable_sites(err_file)
    total_sites <- 10000  # Real value
    
    # Calculate mutation rates per site
    true_mut_rate_per_site <- ifelse(!is.na(true_values$tumourRate), true_values$tumourRate / total_sites, NA)
    est_mut_rate_per_site <- ifelse(!is.na(result_values$tumourRate), result_values$tumourRate / variable_sites, NA)
    
    data.frame(
      model_folder  = model_folder,
      script        = model,
      r_Value       = r_value,
      minD          = params$minD,
      minVar        = params$minVar,
      variable_sites= variable_sites,
      true_tumourAge= true_values$tumourAge,
      true_lambda   = true_values$lambda,
      true_rootRate = true_values$rootRate,
      true_tumourRate= true_values$tumourRate,
      true_mut_rate_per_site = true_mut_rate_per_site,
      est_treeHeight= result_values$tree_height,
      est_lambda    = result_values$lambda,
      est_rootRate  = result_values$rootRate,
      est_tumourRate= result_values$tumourRate,
      est_mut_rate_per_site = est_mut_rate_per_site,
      stringsAsFactors = FALSE
    )
  } else {
    NULL
  }
}

# Main processing loop
all_results <- list()

for (model_folder in model_folders) {
  cat("\n===========================================\n")
  cat("Processing model folder:", model_folder, "\n")
  cat("===========================================\n")
  
  for (model in model_types) {
    cat("  Processing script:", model, "\n")
    
    # Use parallel processing for r values
    n_cores <- detectCores(logical = FALSE)
    model_results <- mclapply(0:99, function(r) {
      process_r_model(r, model, model_folder)
    }, mc.cores = n_cores)
    
    model_results <- Filter(Negate(is.null), model_results)
    all_results   <- c(all_results, model_results)
    
    cat("    Found", length(model_results), "valid results\n")
  }
}

# Combine all results
param_results <- rbindlist(all_results, fill = TRUE)

# Calculate errors/differences from true values
param_results <- param_results %>%
  mutate(
    error_lambda = est_lambda - true_lambda,
    error_rootRate = est_rootRate - true_rootRate,
    error_tumourRate = est_tumourRate - true_tumourRate,
    error_mut_rate_per_site = est_mut_rate_per_site - true_mut_rate_per_site,
    error_abs_lambda = abs(error_lambda),
    error_abs_rootRate = abs(error_rootRate),
    error_abs_tumourRate = abs(error_tumourRate),
    error_abs_mut_rate_per_site = abs(error_mut_rate_per_site),
    pct_error_lambda = ifelse(true_lambda != 0, (error_lambda / abs(true_lambda)) * 100, NA),
    pct_error_rootRate = ifelse(true_rootRate != 0, (error_rootRate / abs(true_rootRate)) * 100, NA),
    pct_error_tumourRate = ifelse(true_tumourRate != 0, (error_tumourRate / abs(true_tumourRate)) * 100, NA),
    pct_error_mut_rate_per_site = ifelse(true_mut_rate_per_site != 0, (error_mut_rate_per_site / abs(true_mut_rate_per_site)) * 100, NA)
  )

# Save the complete results
write.csv(param_results, file = file.path(output_dir, "parameter_analysis_results.csv"), row.names = FALSE)
cat("\n✓ Complete parameter results saved to:", file.path(output_dir, "parameter_analysis_results.csv"), "\n")

# Create summary statistics by model and script
summary_stats <- param_results %>%
  group_by(model_folder, script, minD, minVar) %>%
  summarise(
    n = n(),
    mean_variable_sites = mean(variable_sites, na.rm = TRUE),
    
    mean_est_lambda = mean(est_lambda, na.rm = TRUE),
    sd_est_lambda = sd(est_lambda, na.rm = TRUE),
    mean_true_lambda = mean(true_lambda, na.rm = TRUE),
    mean_error_lambda = mean(error_lambda, na.rm = TRUE),
    sd_error_lambda = sd(error_lambda, na.rm = TRUE),
    mean_pct_error_lambda = mean(pct_error_lambda, na.rm = TRUE),
    
    mean_est_rootRate = mean(est_rootRate, na.rm = TRUE),
    sd_est_rootRate = sd(est_rootRate, na.rm = TRUE),
    mean_true_rootRate = mean(true_rootRate, na.rm = TRUE),
    mean_error_rootRate = mean(error_rootRate, na.rm = TRUE),
    sd_error_rootRate = sd(error_rootRate, na.rm = TRUE),
    mean_pct_error_rootRate = mean(pct_error_rootRate, na.rm = TRUE),
    
    mean_est_tumourRate = mean(est_tumourRate, na.rm = TRUE),
    sd_est_tumourRate = sd(est_tumourRate, na.rm = TRUE),
    mean_true_tumourRate = mean(true_tumourRate, na.rm = TRUE),
    mean_error_tumourRate = mean(error_tumourRate, na.rm = TRUE),
    sd_error_tumourRate = sd(error_tumourRate, na.rm = TRUE),
    mean_pct_error_tumourRate = mean(pct_error_tumourRate, na.rm = TRUE),
    
    mean_est_mut_rate_per_site = mean(est_mut_rate_per_site, na.rm = TRUE),
    sd_est_mut_rate_per_site = sd(est_mut_rate_per_site, na.rm = TRUE),
    mean_true_mut_rate_per_site = mean(true_mut_rate_per_site, na.rm = TRUE),
    mean_error_mut_rate_per_site = mean(error_mut_rate_per_site, na.rm = TRUE),
    sd_error_mut_rate_per_site = sd(error_mut_rate_per_site, na.rm = TRUE),
    mean_pct_error_mut_rate_per_site = mean(pct_error_mut_rate_per_site, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(model_folder, script, minD, minVar)

write.csv(summary_stats, file = file.path(output_dir, "parameter_summary_statistics.csv"), row.names = FALSE)
cat("✓ Summary statistics saved to:", file.path(output_dir, "parameter_summary_statistics.csv"), "\n")

# Print summary to console
cat("\n===========================================\n")
cat("SUMMARY OF PROCESSED DATA\n")
cat("===========================================\n")
cat("Total records:", nrow(param_results), "\n")
cat("Model folders:", paste(unique(param_results$model_folder), collapse = ", "), "\n")
cat("Scripts:", paste(unique(param_results$script), collapse = ", "), "\n")
cat("R values range: 0 -", max(param_results$r_Value, na.rm = TRUE), "\n\n")

cat("Sample of parameter ranges:\n")
cat("Lambda (estimated): ", min(param_results$est_lambda, na.rm=TRUE), " to ", max(param_results$est_lambda, na.rm=TRUE), "\n")
cat("Root Rate (estimated): ", min(param_results$est_rootRate, na.rm=TRUE), " to ", max(param_results$est_rootRate, na.rm=TRUE), "\n")
cat("Tumour Rate (estimated): ", min(param_results$est_tumourRate, na.rm=TRUE), " to ", max(param_results$est_tumourRate, na.rm=TRUE), "\n")
cat("Variable Sites (inferred): ", min(param_results$variable_sites, na.rm=TRUE), " to ", max(param_results$variable_sites, na.rm=TRUE), "\n")
cat("Mutation Rate per Site (estimated): ", min(param_results$est_mut_rate_per_site, na.rm=TRUE), " to ", max(param_results$est_mut_rate_per_site, na.rm=TRUE), "\n")
cat("Mutation Rate per Site (true): ", min(param_results$true_mut_rate_per_site, na.rm=TRUE), " to ", max(param_results$true_mut_rate_per_site, na.rm=TRUE), "\n")

cat("\n✓ Parameter extraction and analysis complete!\n")