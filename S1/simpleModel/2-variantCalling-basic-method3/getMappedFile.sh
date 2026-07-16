#!/bin/bash

# ============================================================
# 9-vcf2fasta.sh
# Run from: 2-variantCalling/
#
# For each DxVx/merge/root_ref_<suffix>.fasta:
#   - Find matching DxVx/merge/variants_<suffix>.vcf
#   - Run vcf2fasta.py
# Output: DxVx/merge/raw_alignment_<suffix>.fasta
# ============================================================

BASE_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "================================================"
echo "Working directory: ${BASE_DIR}"
echo "Start time: $(date)"
echo "================================================"

convert_count=0
warn_count=0

for dxvx_dir in "${BASE_DIR}"/D*V*/; do
    [[ -d "${dxvx_dir}" ]] || continue
    sample=$(basename "${dxvx_dir}")
    
    merge_dir="${dxvx_dir}merge"
    output_dir="${dxvx_dir}alignment"
    mkdir -p ${output_dir}
    
    if [[ ! -d "${merge_dir}" ]]; then
        echo " [SKIP] ${sample}: no merge/ directory"
        continue
    fi

    echo ""
    echo "Processing: ${sample}"

    for ref_seq in "${merge_dir}"/root_ref_*.fasta; do
        [[ -f "${ref_seq}" ]] || continue

        # root_ref_r0_d1_v1.fasta -> r0_d1_v1
        suffix=$(basename "${ref_seq}" .fasta | sed 's/root_ref_//')

        vcf_file="${merge_dir}/variants_${suffix}.vcf"
        output="${output_dir}/alignment_${suffix}.fasta"

        if [[ -f "${vcf_file}" ]]; then
            echo "  ${suffix}: ${vcf_file} -> alignment_${suffix}.fasta"
            python ~/vcf2fasta/src/vcf2fasta.py \
                "${vcf_file}" "${output}" \
                --ref "${ref_seq}" \
                --refType nd16
            if [[ $? -eq 0 ]]; then
                echo " OK"
                convert_count=$((convert_count + 1))
            else
                echo " ERROR running vcf2fasta for ${suffix}" >&2
                warn_count=$((warn_count + 1))
            fi
        else
            echo " [WARN] VCF not found for ${suffix}: ${vcf_file}"
            warn_count=$((warn_count + 1))
        fi
    done
done

echo ""
echo "================================================"
echo "Done. End time: $(date)"
echo "  Converted : ${convert_count}"
echo "  Warnings  : ${warn_count}"
echo "================================================"
