#!/bin/bash

# Input and output directories
script_folder="../1-simulation/pileup"
output_folder="./filtered"

# Recreate output directory
rm -rf "$output_folder"
mkdir -p "$output_folder"

# Loop through all subdirectories
for subfolder in "${script_folder}"/*/; do
    [ -d "$subfolder" ] || continue

    subfolder_name=$(basename "$subfolder")
    echo "Processing subfolder: $subfolder_name"

    # Create matching output subdirectory
    out_subfolder="${output_folder}/${subfolder_name}"
    mkdir -p "$out_subfolder"

    # Process all VCF files
    for file in "${subfolder}"/*.vcf; do
        [ -f "$file" ] || continue

        name=$(basename "$file")
        output="${out_subfolder}/${name}"

        echo "Creating filtered output file: $output"

        # Write header lines first, then filtered variant lines
        awk '
        /^#/ {print; next}
        $7 == "." && length($5) == 1 {print}
        ' "$file" > "$output"

    done
done