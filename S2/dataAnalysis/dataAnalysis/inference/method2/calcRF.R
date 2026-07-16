# install these packages first
# Load necessary packages
library(ape)
library(phangorn)

# Define the model folder names
model_folders <- c("simpleModel", "complexModel", "calibratedModel")

# Initialize an empty list to store all results
all_results <- list()

# Define the base paths
base_path <- "../../../"

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

# Loop through each model folder
for (model_folder in model_folders) {
  
  # Define folder paths for the current model
  trueTreeFolder <- file.path(base_path, model_folder, "1-simulation/trees")
  resultFolder <- file.path(base_path, model_folder, "3-xml-method2/results")
  paramsFolder <- file.path(base_path, model_folder, "1-simulation/pileup")
  
  cat("Processing model:", model_folder, "\n")
  cat("True tree folder:", trueTreeFolder, "\n")
  cat("Result folder:", resultFolder, "\n")
  cat("(Full analysis - all tips included)\n\n")
  
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
        
        # Ensure that the tip labels are the same
        if (!all(sort(true_tree$tip.label) == sort(consensus_tree$tip.label))) {
          cat("Warning: Labels do not match between true tree and consensus tree for", 
              true_tree_name, "with model", model, "\n")
          
          # Reorder the consensus tree's tip labels to match true tree
          consensus_tree <- ape::keep.tip(consensus_tree, true_tree$tip.label)
        }
        
        # Calculate the RF distance between the true tree and the consensus tree
        rf_distance <- phangorn::RF.dist(tree1 = true_tree, tree2 = consensus_tree, 
                                          normalize = TRUE, check.labels = TRUE, rooted = TRUE)
        
        # Calculate the tree height of consensus and true trees
        tree_height_true <- max(ape::node.depth.edgelength(true_tree))
        tree_height_consensus <- max(ape::node.depth.edgelength(consensus_tree))
        
        # Calculate the total tree length (sum of branch lengths)
        tree_length_consensus <- sum(consensus_tree$edge.length)
        tree_length_true <- sum(true_tree$edge.length)
        
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
          "num_tips" = length(true_tree$tip.label),
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
  output_file <- file.path(output_folder, "tree_analysis_results_full_analysis.csv")
  
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
  
  # Optionally, print a message confirming the file has been saved
  cat("\nResults saved to:", output_file, "\n")
}