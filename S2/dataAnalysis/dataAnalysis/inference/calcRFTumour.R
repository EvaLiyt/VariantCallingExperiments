# install these packages first
# Load necessary packages
library(ape)
library(phangorn)
library(dplyr)

# Define the model folder names
model_folders <- c("simpleModel", "complexModel", "calibratedModel")

# Initialize an empty list to store all results
all_results <- list()

# Define the base paths
base_path <- "../../"

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

# Define the consensus model types
model_types <- c("ErrModel", "wtErrModel", "TrueErrModel", "TruewtErrModel")

# Function to filter tree to only tumour clades (clade_ prefix)
# Limits to maximum 40 clades
filter_tree_to_40_tumour_clades <- function(tree, model_folder) {
  
  # For simpleModel, keep all tips
  if (model_folder == "simpleModel") {
    return(tree)
  }
  
  # For other models, keep only tips with clade_ prefix
  clade_tips <- tree$tip.label[grepl("^clade_", tree$tip.label)]
  
  if (length(clade_tips) == 0) {
    cat("Warning: No clade_ prefixed tips found in tree\n")
    return(NULL)
  }
  
  # Limit to maximum 40 clades
  if (length(clade_tips) > 40) {
    cat("  Note: Found", length(clade_tips), "clades, using first 40\n")
    clade_tips <- clade_tips[1:40]
  }
  
  # Keep only the clade tips
  filtered_tree <- ape::keep.tip(tree, clade_tips)
  
  return(filtered_tree)
}

# Loop through each model folder
for (model_folder in model_folders) {
  
  # Define folder paths for the current model
  trueTreeFolder <- file.path(base_path, model_folder, "1-simulation/trees")
  resultFolder <- file.path(base_path, model_folder, "3-xml/results")
  paramsFolder <- file.path(base_path, model_folder, "1-simulation/pileup")
  
  cat("Processing model:", model_folder, "\n")
  cat("True tree folder:", trueTreeFolder, "\n")
  cat("Result folder:", resultFolder, "\n")
  if (model_folder != "simpleModel") {
    cat("(Filtering to maximum 40 tumour clades with 'clade_' prefix)\n")
  } else {
    cat("(Keeping all tips for simpleModel)\n")
  }
  cat("\n")
  
  # Check if the required folders exist
  if (!dir.exists(trueTreeFolder)) {
    cat("Warning: True tree folder does not exist:", trueTreeFolder, "\n\n")
    next
  }
  
  if (!dir.exists(resultFolder)) {
    cat("Warning: Result folder does not exist:", resultFolder, "\n\n")
    next
  }
  
  # List all true tree files (specifically _psi.trees files)
  true_tree_files <- list.files(trueTreeFolder, pattern = "_psi\\.trees$", full.names = TRUE)
  if(model_folder=="calibratedModel"){
    true_tree_files <- list.files(trueTreeFolder, pattern = "_tree\\.trees$", full.names = TRUE)
  }
  
  
  if (length(true_tree_files) == 0) {
    cat("No tree files found in:", trueTreeFolder, "\n\n")
    next
  }
  
  cat("Found", length(true_tree_files), "tree files\n\n")
  
  # Loop through each true tree file
  for (true_tree_file in true_tree_files) {
    # Extract the base name of the true tree file
    true_tree_name <- tools::file_path_sans_ext(basename(true_tree_file))
    
    # Extract the 'r' value (e.g., r0 from sequenceGeneration_r0_psilabelled.trees)
    r_value <- sub(".*_r([0-9]+).*", "\\1", true_tree_name)
    
    # Read the true tree
    true_tree <- tryCatch(
      ape::read.nexus(true_tree_file),
      error = function(e) {
        cat("Error reading tree file:", true_tree_file, "\n")
        return(NULL)
      }
    )
    
    if (is.null(true_tree)) next
    
    # Filter true tree to 40 tumour clades
    true_tree_filtered <- filter_tree_to_40_tumour_clades(true_tree, model_folder)
    if (is.null(true_tree_filtered)) next
    
    # Loop through each consensus model type
    for (model in model_types) {
      consensus_tree_file <- file.path(resultFolder, paste0(model, "_r", r_value, ".trees.txt"))
      
      # Check if the consensus tree file exists
      if (file.exists(consensus_tree_file)) {
        # Read the consensus tree
        consensus_tree <- tryCatch(
          ape::read.nexus(consensus_tree_file),
          error = function(e) {
            cat("Error reading consensus tree file:", consensus_tree_file, "\n")
            return(NULL)
          }
        )
        
        if (is.null(consensus_tree)) next
        
        # Filter consensus tree to 40 tumour clades
        consensus_tree_filtered <- filter_tree_to_40_tumour_clades(consensus_tree, model_folder)
        if (is.null(consensus_tree_filtered)) next
        
        # Ensure that the tip labels are the same between filtered trees
        if (!all(sort(true_tree_filtered$tip.label) == sort(consensus_tree_filtered$tip.label))) {
          cat("Warning: Labels do not match between true tree and consensus tree for", 
              true_tree_name, "with model", model, "\n")
          
          # Reorder the consensus tree's tip labels to match true tree
          consensus_tree_filtered <- ape::keep.tip(consensus_tree_filtered, true_tree_filtered$tip.label)
        }
        
        # Calculate the normalized RF distance between the filtered trees
        rf_distance <- phangorn::RF.dist(tree1 = true_tree_filtered, tree2 = consensus_tree_filtered, 
                                          normalize = TRUE, check.labels = TRUE, rooted = TRUE)
        
        # Calculate the tree height of consensus and true trees (filtered)
        tree_height_true <- max(ape::node.depth.edgelength(true_tree_filtered))
        tree_height_consensus <- max(ape::node.depth.edgelength(consensus_tree_filtered))
        
        # Calculate the total tree length (sum of branch lengths)
        tree_length_consensus <- sum(consensus_tree_filtered$edge.length)
        tree_length_true <- sum(true_tree_filtered$edge.length)
        
        # Mean Absolute Error (MAE) for tree height
        mae_height <- abs(tree_height_true - tree_height_consensus)
        
        # Mean Absolute Error (MAE) for tree length
        mae_length <- abs(tree_length_true - tree_length_consensus)
        
        # Parse minD and minVar
        param_file <- file.path(paramsFolder, paste0("r", r_value), paste0("parameters_r", r_value, ".txt"))
        params <- parse_params(param_file)
        
        # Access minD and minVar
        minD <- params$minD
        minVar <- params$minVar
        
        # Store the results in a data frame
        result <- data.frame(
          "model_folder" = model_folder,
          "true_tree" = true_tree_name,
          "script" = model,
          "r_Value" = r_value,
          "minD" = minD,
          "minVar" = minVar,
          "num_clades" = length(true_tree_filtered$tip.label),
          "RF_Dist" = rf_distance,
          "true_treeH" = tree_height_true,
          "consensus_treeH" = tree_height_consensus,
          "true_treeL" = tree_length_true,
          "consensus_treeL" = tree_length_consensus,
          "MAE_treeH" = mae_height,
          "MAE_treeL" = mae_length,
          stringsAsFactors = FALSE
        )
        
        # Append the result to the list of all results
        all_results[[length(all_results) + 1]] <- result
        
      } else {
        # If the file doesn't exist, print a message for debugging
        cat("File not found:", consensus_tree_file, "\n")
      }
    }
  }
  
  cat("\n")
}

# Check if any results are found
if (length(all_results) == 0) {
  cat("No results found.\n")
} else {
  # Combine all results into one data frame
  final_results <- do.call(rbind, all_results)
  
  # Reset row names for cleaner output
  rownames(final_results) <- NULL
  
  # Define the desired folder path for saving the results
  output_folder <- "./summary"
  
  # Ensure the folder exists, if not, create it
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  
  # Define the full path for the output file
  output_file <- file.path(output_folder, "tree_analysis_results_40_clades.csv")
  
  # Save the final results as a CSV file
  write.csv(final_results, output_file, row.names = FALSE)
  
  # Optionally, print a summary of the results
  cat("Results summary:\n")
  cat("Total rows processed:", nrow(final_results), "\n")
  cat("Models compared:", paste(unique(final_results$model_folder), collapse = ", "), "\n")
  cat("Consensus models:", paste(unique(final_results$script), collapse = ", "), "\n")
  cat("\n")
  
  # Print first few rows
  cat("First few rows of results:\n")
  print(head(final_results, 10))
  
  # Print number of clades by model
  cat("\nNumber of tumour clades by model:\n")
  clade_counts <- final_results %>%
    distinct(model_folder, num_clades) %>%
    arrange(model_folder)
  print(clade_counts)
  
  # Optionally, print a message confirming the file has been saved
  cat("\nResults saved to:", output_file, "\n")
}