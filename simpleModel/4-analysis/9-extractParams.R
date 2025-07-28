library(data.table)
library(parallel)

# Define the folder paths
trueLogFolder <- "./trueLogs"
resultFolder <- "./parsed_results"
paramsFolder <- "../1-simulation/pileup"

# Create a list of model types
model_types <- c("ErrModel", "wtErrModel", "TrueErrModel", "TruewtErrModel")

# Ensure output directory exists
if (!dir.exists("./summary")) dir.create("./summary", recursive = TRUE)

# Efficient helper to extract values from the true log files
extract_true_values <- function(filepath) {
  lines <- readLines(filepath, warn = FALSE)
  
  # Remove blank and comment lines
  lines <- lines[grepl("\\S", lines) & !grepl("^#", lines)]
  if (length(lines) < 2) {
    warning(paste("File format error:", filepath))
    return(list(tumourAge=NA, lambda=NA, rootRate=NA, tumourRate=NA, rootAge=NA))
  }
  
  # Detect delimiter
  if (grepl("\t", lines[1])) {
    delimiter <- "\t"
  } else {
    delimiter <- "[[:space:]]+"
  }
  
  # Split header and value lines
  params <- unlist(strsplit(trimws(lines[1]), delimiter))
  values <- as.numeric(unlist(strsplit(trimws(lines[2]), delimiter)))
  
  if (length(params) != length(values)) {
    warning(paste("Parameter/value count mismatch in", filepath))
    return(list(tumourAge=NA, lambda=NA, rootRate=NA, tumourRate=NA, rootAge=NA))
  }
  
  get_param_value <- function(param_prefix) {
    matches <- grep(param_prefix, params, fixed = TRUE)
    if (length(matches) > 0) {
      return(values[matches[1]])
    } else {
      return(NA)
    }
  }
  
  lambda <- get_param_value("lambda")
  tumourAge <- get_param_value("tumourAge")
  rootAge <- get_param_value("rootAge")
  rootRate <- get_param_value("rootRate")
  
  # Flexible tumourRate extraction
  tumourRate <- get_param_value("tumourRate")
  if (is.na(tumourRate)) {
    tumourRate <- get_param_value("mu")
  }
  
  return(list(
    tumourAge = tumourAge,
    lambda = lambda,
    rootRate = rootRate,
    tumourRate = tumourRate,
    rootAge = rootAge
  ))
}

# Helper function to extract values from result files
extract_result_values <- function(filepath) {
  result_data <- read.table(filepath, header = TRUE, stringsAsFactors = FALSE)
  params <- colnames(result_data)
  means <- as.numeric(result_data[1, ])
  
  get_param_value <- function(param_name) {
    if (param_name %in% params) {
      return(means[params == param_name])
    } else {
      return(NA)
    }
  }
  
  lambda <- get_param_value("lambda")
  rootRate <- get_param_value("rootRate")
  tree_height <- get_param_value("tree.height")
  tumourRate <- get_param_value("tumourRate")
  if (is.na(tumourRate)) tumourRate <- get_param_value("mu")
  
  return(list(lambda = lambda, rootRate = rootRate, tumourRate = tumourRate, tree_height = tree_height))
}

# Function to process a single r/model combination
process_r_model <- function(r_value, model) {
  true_log_file <- file.path(trueLogFolder, paste0("sequenceGeneration_r", r_value, ".log"))
  result_file <- file.path(resultFolder, paste0(model, "_r", r_value, "_stats.log"))
  
  if (file.exists(true_log_file) && file.exists(result_file)) {
    true_values <- extract_true_values(true_log_file)
    result_values <- extract_result_values(result_file)
    
    return(data.frame(
      script = model,
      r_Value = r_value,
      true_tumourAge = true_values$tumourAge,
      true_lambda = true_values$lambda,
      true_rootRate = true_values$rootRate,
      true_tumourRate = true_values$tumourRate,
      est_treeHeight = result_values$tree_height,
      est_lambda = result_values$lambda,
      est_rootRate = result_values$rootRate,
      est_tumourRate = result_values$tumourRate,
      stringsAsFactors = FALSE
    ))
  } else {
    return(NULL)
  }
}

# Detect number of cores
n_cores <- detectCores(logical = FALSE)

# Collect results
all_results <- list()

# Process each model in parallel
for (model in model_types) {
  cat("Processing model:", model, "\n")
  
  model_results <- mclapply(0:99, function(r) process_r_model(r, model), mc.cores = n_cores)
  model_results <- Filter(Negate(is.null), model_results)
  
  all_results <- c(all_results, model_results)
}

# Combine all results
final_results <- rbindlist(all_results, fill = TRUE)

# Save results
write.csv(final_results, file = "./summary/parameter_results.csv", row.names = FALSE)

cat("Final results (first 10 rows):\n")
print(head(final_results, 10))