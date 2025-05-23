#!/bin/bash

# Load modules
module load GCCcore/11.3.0
module load BCFtools/1.19-GCC-11.3.0

# Define input/output directories
out_dir="./merge"
vcf_dir="./changedVCF"

# Reset output directory
rm -rf "$out_dir"
mkdir -p "$out_dir"

# Traverse all VCF files and group them by 'r'
for vcf_file in "$vcf_dir"/*/*.vcf.gz.gz; do
    if [[ -f "$vcf_file" ]]; then
        # Extract r value from file name (e.g., r0 from r0_0.vcf.gz.gz)
        file_name=$(basename "$vcf_file")
        r_value=$(echo "$file_name" | cut -d'_' -f1)  # r0

        # Build list file for this group
        vcf_list="${out_dir}/vcf_list_${r_value}.txt"
        echo "$vcf_file" >> "$vcf_list"
    fi
done

# Merge files in each list
for vcf_list in "$out_dir"/vcf_list_r*.txt; do
    r_value=$(basename "$vcf_list" .txt | sed 's/vcf_list_//')
    out_vcf="${out_dir}/merged_${r_value}.vcf"

    # Merge using bcftools
    bcftools merge --file-list "$vcf_list" -o "$out_vcf"

    if [ $? -eq 0 ]; then
        echo "Successfully merged VCFs for ${r_value} into ${out_vcf}"
    else
        echo "Error merging for ${r_value}"
        exit 1
    fi
done

echo "Done merging all VCF files by r group"
