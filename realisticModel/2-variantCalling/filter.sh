#!/bin/bash

# Define the input folder containing subfolders and the output folder
script_folder="../1-simulation/pileup"
output_folder="./filtered"

# Remove the output folder if it exists, then recreate it
rm -rf "$output_folder"
mkdir -p "$output_folder"

# Iterate over all subfolders in the input folder
for subfolder in "${script_folder}"/*/; do
    # Ensure that it is indeed a directory
    if [ -d "$subfolder" ]; then
        # Get the name of the subfolder (e.g., r0, r1, ..., r99)
        subfolder_name=$(basename "$subfolder")
        echo "Processing subfolder: $subfolder_name"
        
        # Create a corresponding output subfolder for each input subfolder
        mkdir -p "${output_folder}/${subfolder_name}"

        # Loop over each VCF file in the current subfolder
        for file in "${subfolder}"/*.vcf; do
            # Only proceed if the file exists (in case there are no VCF files)
            if [ -f "$file" ]; then
                # Strip 'sorted_' prefix and '.bam.vcf' suffix from the filename
                name="${file##*/}"        # Get the filename without the path

                # Define the output file path
                output="${output_folder}/${subfolder_name}/${name}"

                # Print status message for debugging
                echo "Creating filtered output file: $output"

                # Filter lines with '#' or with missing data '.' in the 7th column and write to the output file
                awk -F'\t' '/^#/ || ($7 == "." && $5 !~ /[\+\-]/) { print }' "$file" > "$output"


            fi
        done
    fi
done
