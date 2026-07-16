#!/usr/bin/env bash

set -euo pipefail

# Input folder containing replicate subfolders
script_folder="step3"

# Output folder
output_folder="./filtered"

# Remove and recreate output directory
rm -rf "$output_folder"
mkdir -p "$output_folder"

# Iterate through all subfolders (e.g. r0, r1, ...)
for subfolder in "${script_folder}"/*/; do
    # Skip if not a directory
    [[ -d "$subfolder" ]] || continue

    subfolder_name=$(basename "$subfolder")
    echo "Processing subfolder: ${subfolder_name}"

    # Create matching output subfolder
    mkdir -p "${output_folder}/${subfolder_name}"

    # Process each VCF file
    for file in "${subfolder}"/*.vcf; do
        # Skip if no VCF files exist
        [[ -f "$file" ]] || continue

        name=$(basename "$file")
        output="${output_folder}/${subfolder_name}/${name}"

        echo "  Filtering: ${name}"

        # Keep:
        #   - header lines
        #   - records with FILTER="."
        #   - single-base ALT alleles only (no indels)
        awk '
            /^#/ {
                print
                next
            }

            $7 == "." && length($5) == 1 {
                print
            }
        ' "$file" > "$output"
    done
done

echo "Filtering complete."