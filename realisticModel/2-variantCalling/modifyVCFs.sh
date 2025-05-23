#!/bin/bash

# Define the input and output directories on NeSI
input_folder="./filtered"
output_folder="./modifiedVCF"

# Remove and recreate output folder
rm -rf "$output_folder"
mkdir -p "$output_folder"

# Traverse all subfolders (r0, r1, ..., r99) in the input folder
for subfolder in "${input_folder}"/*/; do
    # Ensure that it is indeed a directory
    if [ -d "$subfolder" ]; then
        # Get the name of the subfolder (e.g., r0, r1, ..., r99)
        subfolder_name=$(basename "$subfolder")
        echo "Processing subfolder: $subfolder_name"
        
        # Create the corresponding output subfolder
        mkdir -p "${output_folder}/${subfolder_name}"
        
        # Traverse all VCF files in the subfolder
        for input_file in "$subfolder"/*.vcf; do
            if [ -f "$input_file" ]; then
                # Get the filename and define the output file path
                filename=$(basename "$input_file")
                output_file="${output_folder}/${subfolder_name}/$filename"
                
                # Modify the VCF file
                # 1. Delete lines starting with "##FILTER=<ID=No.variants<2"
                # 2. Replace "Number=G" with "Number=1"
                sed -e '/^##FILTER=<ID=No\.variants<2/d' \
                    -e 's/Number=G/Number=1/' "$input_file" > "$output_file"


                echo "Processed $filename"
            fi
        done
    fi
done

echo "Finished modifying all VCF files."
