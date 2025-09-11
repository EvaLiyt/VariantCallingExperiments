# install these packages first
# Load necessary packages
library(ape)
library(phangorn)

# Define the folder paths
trueTreeFolder <- "../1-simulation/trees"
resultFolder <- "../3-xml/results"
paramsFolder <- "../1-simulation/pileup"

# Create a list of minD values and model types
model_types <- c("ErrModel", "wtErrModel", "TrueErrModel", "TruewtErrModel")
# Initialize an empty list to store all results
all_results <- list()

# List all true tree files (assuming they are in .trees format)
true_tree_files <- list.files(trueTreeFolder, pattern = "\\.trees$", full.names = TRUE)

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

# Loop through each true tree file
for (true_tree_file in true_tree_files) {
  # Extract the base name of the true tree file (e.g., sequenceGeneration_r0_psilabelled.trees)
  true_tree_name <- tools::file_path_sans_ext(basename(true_tree_file))
  
  # Extract the 'r' value (e.g., r0 from sequenceGeneration_r0_psilabelled.trees)
  r_value <- sub(".*_r([0-9]+).*", "\\1", true_tree_name)
  
  # Read the true tree
  true_tree <- ape::read.nexus(true_tree_file)
  
  # Loop through each combination of model and minD
  for (model in model_types) {
    consensus_tree_file <- file.path(resultFolder, paste0(model, "_r", r_value,".trees.txt"))
      
    # Check if the consensus tree file exists
    if (file.exists(consensus_tree_file)) {
      # Read the consensus tree
      consensus_tree <- ape::read.nexus(consensus_tree_file)
      
      # Ensure that the tip labels are the same
      if (!all(sort(true_tree$tip.label) == sort(consensus_tree$tip.label))) {
        cat("Warning: Labels do not match between true tree and consensus tree for", true_tree_name, "with model", model, "and minD", minD, "\n")
        
        # Optionally, you can reorder the trees' tip labels to match
        # Reordering the consensus tree's tip labels to match true tree
        consensus_tree <- keep.tip(consensus_tree, true_tree$tip.label)
      }
        
        # Calculate the RF distance between the true tree and the consensus tree
        rf_distance <- RF.dist(tree1 = true_tree, tree2 = consensus_tree, normalize = TRUE, check.labels = TRUE, rooted = TRUE)
        
        # Calculate the tree height of consensus and true trees
        tree_height_true <- max(node.depth.edgelength(true_tree))
        tree_height_consensus <- max(node.depth.edgelength(consensus_tree))
        
        # Calculate the total tree length (sum of branch lengths)
        tree_length_consensus <- sum(consensus_tree$edge.length)
        tree_length_true <- sum(true_tree$edge.length)
        
        # Mean Absolute Error (MAE) for tree height
        mae_height <- abs(tree_height_true - tree_height_consensus)
        
        # Mean Absolute Error (MAE) for tree length
        mae_length <- abs(tree_length_true - tree_length_consensus)
        
        # minD and minVar
        param_file <- file.path(paramsFolder, paste0("r", r_value), paste0("parameters_r", r_value, ".txt"))
        params <- parse_params(param_file)
        
        # Access minD and minVar
        minD <- params$minD
        minVar <- params$minVar
        
        # Store the results in a data frame
        result <- data.frame(
          "true_tree" = true_tree_name,
          "script" = model,
          "r_Value" = r_value,
          "minD" = minD,
          "minVar" = minVar,
          "RF_Dist" = rf_distance,
          "true_treeH" = tree_height_true,
          "consensus_treeH" = tree_height_consensus,
          "true_treeL" = tree_length_true,
          "consensus_treeL" = tree_length_consensus,
          "MAE_treeH" = mae_height,
          "MAE_treeL" = mae_length
        )
        
        # Append the result to the list of all results
        all_results[[length(all_results) + 1]] <- result
      } else {
        # If the file doesn't exist, print a message for debugging
        cat("File not found:", consensus_tree_file, "\n")
      }
    }
  }

# Check if any results are found
if (length(all_results) == 0) {
  cat("No results found.\n")
} else {
  # Combine all results into one data frame
  final_results <- do.call(rbind, all_results)
  
  # Define the desired folder path for saving the results
  output_folder <- "./summary"
  
  # Ensure the folder exists, if not, create it
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  
  # Define the full path for the output file
  output_file <- file.path(output_folder, "tree_analysis_results.csv")
  
  # Save the final results as a CSV file
  write.csv(final_results, output_file, row.names = FALSE)
  
  # Optionally, print a message confirming the file has been saved
  cat("Results saved to:", output_file, "\n")
  
}
