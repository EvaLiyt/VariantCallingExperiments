#!/bin/bash
#SBATCH -J RM_vcf2fasta
#SBATCH --time=5:00:00        # Adjust walltime as needed

# Set paths
file_folder="./merge"

# Loop through each reference fasta file (new_ref_rX.fasta)
for ref_seq in ${file_folder}/root_ref_r*.fasta; do
    # Extract the rX part from the file name
    suffix=$(basename "$ref_seq" .fasta | sed 's/root_ref_//')

    # Construct the corresponding VCF file path (e.g., variants_rX.vcf)
    vcf_file="${file_folder}/combined_${suffix}.vcf"
    
    # Check if the corresponding VCF file exists
    if [ -f "$vcf_file" ]; then
        # Construct the output FASTA file path
        output="${file_folder}/raw_alignment_${suffix}.fasta"

        # Run the vcf2fasta tool
        python ~/vcf2fasta/src/vcf2fasta.py "$vcf_file" "$output" --ref "$ref_seq" --refType nd16

        echo "VCF to FASTA conversion completed for ${vcf_file} -> ${output}"
    else
        echo "Warning: VCF file ${vcf_file} not found for reference ${ref_seq}"
    fi
done

