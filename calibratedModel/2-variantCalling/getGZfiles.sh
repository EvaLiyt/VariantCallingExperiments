#!/bin/bash

# Load the bcftools module
module load GCCcore/11.3.0
module load BCFtools/1.19-GCC-11.3.0

# Define the input and reference sequence
input_folder="./modifiedVCF"

# Prepare for processing
for subfolder in "${input_folder}"/*; do
    # Ensure that it is indeed a directory
    if [ -d "$subfolder" ]; then
        # Get the name of the subfolder (e.g., r0, r1, ..., r99)
        subfolder_name=$(basename "$subfolder")
        echo "Processing subfolder: $subfolder_name"
        
        for file in "${subfolder}"/*vcf; do
            filename=$(basename "${file}")
            
            # Convert to gz file
            gz_file="${subfolder}/${filename%.vcf}.vcf.gz"  # Ensure proper naming
            bcftools view "${file}" -Oz -o "${gz_file}"

            # Normalisation
            norm_file="${subfolder}/${filename%.vcf}.norm.vcf.gz"  # Avoid double .gz
            bcftools norm -Oz -f ~/data/10K.fna -o "${norm_file}" "${gz_file}"

            # Index the gz file
            bcftools index "${norm_file}"
        done
    fi
done

echo "Done indexing and normalising VCF files."
