#!/usr/bin/env Rscript

library(data.table)
library(parallel)

# -----------------------------
# Paths
# -----------------------------
trueLogFolder <- "../1-simulation/simulationData"
resultFolder  <- "./parsed_results"
paramsFolder  <- "../1-simulation/pileup"
xmlsFolder    <- "../3-xml/xmls"

output_dir <- "./summary"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -----------------------------
# Part 1: Parameter extraction
# -----------------------------
model_types <- c("ErrModel", "wtErrModel", "TrueErrModel", "TruewtErrModel")

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

process_r_model <- function(r_value, model) {
  true_log_file <- file.path(trueLogFolder, paste0("sequenceGeneration_r", r_value, ".log"))
  result_file   <- file.path(resultFolder, paste0(model, "_r", r_value, "_stats.log"))
  if (file.exists(true_log_file) && file.exists(result_file)) {
    true_values   <- extract_true_values(true_log_file)
    result_values <- extract_result_values(result_file)
    data.frame(
      script        = model,
      r_Value       = r_value,
      true_tumourAge= true_values$tumourAge,
      true_lambda   = true_values$lambda,
      true_rootRate = true_values$rootRate,
      true_tumourRate= true_values$tumourRate,
      est_treeHeight= result_values$tree_height,
      est_lambda    = result_values$lambda,
      est_rootRate  = result_values$rootRate,
      est_tumourRate= result_values$tumourRate,
      stringsAsFactors = FALSE
    )
  } else NULL
}

n_cores <- detectCores(logical = FALSE)
all_results <- list()
for (model in model_types) {
  cat("Processing model:", model, "\n")
  model_results <- mclapply(0:99, function(r) process_r_model(r, model), mc.cores = n_cores)
  model_results <- Filter(Negate(is.null), model_results)
  all_results   <- c(all_results, model_results)
}
param_results <- rbindlist(all_results, fill = TRUE)

# -----------------------------
# Part 2: Sites extraction
# -----------------------------
out_files <- list.files(xmlsFolder, pattern="*out$", full.names=TRUE)

sites_results <- lapply(out_files, function(f) {
  filename <- basename(f)
  filename_noext <- sub("\\.out$", "", filename)
  
  parts <- strsplit(filename_noext, "_r")[[1]]
  
  if (length(parts) >= 2) {
    script <- paste(parts[-length(parts)], collapse="_") 
    r_Value <- as.integer(parts[length(parts)])
  } else {
    script <- filename_noext
    r_Value <- NA
  }
  
  lines <- readLines(f, warn=FALSE)
  site_line <- grep("Alignment\\(A\\):", lines, value=TRUE)
  
  if (length(site_line) > 0) {
    sites <- as.integer(sub(".*\\[([0-9]+), ([0-9]+), ([0-9]+)\\].*", "\\3", site_line[1]))
  } else {
    sites <- NA
  }
  
  data.frame(script=script, r_Value=r_Value, sites=sites, stringsAsFactors=FALSE)
})

sites_results <- rbindlist(sites_results, fill=TRUE)


# -----------------------------
# Merge and Save
# -----------------------------
final_results <- merge(param_results, sites_results, by=c("script","r_Value"), all.x=TRUE)

write.csv(final_results, file = file.path(output_dir, "parameter_results.csv"), row.names = FALSE)

cat("Final merged results with sites saved.\n")
