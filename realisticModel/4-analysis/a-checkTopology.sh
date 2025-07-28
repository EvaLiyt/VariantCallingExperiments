#!/bin/bash

# Define folder paths
data_trees_folder="../1-simulation/trees"  # Folder containing the true tree files
results_folder="../3-xml/results"       # Folder containing the MCMC sample tree files
jar_path="/home/yil218/TreeCredibleLevel.jar"  # Path to the TreeCredibleLevel JAR file
output_folder="./topologyCheck"  # Folder for output files

# Set burn-in percentage
burnin=10

# Create the output folder if it doesn't exist
mkdir -p "$output_folder"

# Loop through all .trees files in the results folder
for mcmc_file in "$results_folder"/*.trees; do
    # Extract the filename and r number from the MCMC file
    mcmc_filename=$(basename "$mcmc_file")
    r_number=$(echo "$mcmc_filename" | grep -oP "_r\K[0-9]+") # Extract the number after "_r"

    # Construct the corresponding true tree file path
    true_tree_filename="sequenceGeneration_r${r_number}_psilabelled.trees"
    true_tree_path="$data_trees_folder/$true_tree_filename"

    # Check if the true tree file exists
    if [[ ! -f "$true_tree_path" ]]; then
        echo "True tree file not found for r${r_number}: $true_tree_path"
        continue
    fi

    # Set the output file path, replacing .trees with .txt
    output_file="$output_folder/${mcmc_filename%.trees}.txt"

    # Construct and run the Java command
    echo "Processing: $mcmc_filename"
    echo "Running: java -jar $jar_path -trees $mcmc_file -burnin $burnin -tree $true_tree_path -out $output_file"

    ~/zulu17.54.21-ca-fx-jdk17.0.13-linux_x64/bin/java -jar "$jar_path" \
        -trees "$mcmc_file" \
        -burnin "$burnin" \
        -tree "$true_tree_path" \
        -out "$output_file"

    # Check if the Java command was successful
    if [[ $? -ne 0 ]]; then
        echo "Error processing $mcmc_filename"
    else
        echo "Successfully processed $mcmc_filename"
    fi
done

# Print completion message
echo "Processing complete."
