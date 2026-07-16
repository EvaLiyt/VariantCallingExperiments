#!/bin/bash

# ============================================================
# 6-modifyVCF.sh
# Run from: 2-variantCalling/
#
# For each DxVx/filtered/rx/*.vcf:
#   1. Delete lines starting with "##FILTER=<ID=No.variants"
#   2. Replace "Number=G" with "Number=1"
# Output: DxVx/modifiedVCF/rx/<same filename>
# ============================================================

BASE_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "================================================"
echo "Working directory: ${BASE_DIR}"
echo "Start time: $(date)"
echo "================================================"

total=0
success=0

for dxvx_dir in "${BASE_DIR}"/D*V*/; do
    [[ -d "${dxvx_dir}" ]] || continue
    sample=$(basename "${dxvx_dir}")

    filtered_dir="${dxvx_dir}filtered"
    if [[ ! -d "${filtered_dir}" ]]; then
        echo "  [SKIP] ${sample}: no filtered/ directory"
        continue
    fi

    echo ""
    echo "Processing: ${sample}"

    # Remove and recreate modifiedVCF for this sample
    mod_base="${dxvx_dir}modifiedVCF"
    rm -rf "${mod_base}"
    mkdir -p "${mod_base}"

    for rx_dir in "${filtered_dir}"/r*/; do
        [[ -d "${rx_dir}" ]] || continue
        rep=$(basename "${rx_dir}")

        out_dir="${mod_base}/${rep}"
        mkdir -p "${out_dir}"

        for input_file in "${rx_dir}"/*.vcf; do
            [[ -f "${input_file}" ]] || continue

            fname=$(basename "${input_file}")
            output_file="${out_dir}/${fname}"

            total=$((total + 1))

            sed -e '/^##FILTER=<ID=No\.variants/d' \
                -e 's/Number=G/Number=1/g' \
                "${input_file}" > "${output_file}"

            echo "  ${sample}/modifiedVCF/${rep}/${fname}"
            success=$((success + 1))
        done
    done
done

echo ""
echo "================================================"
echo "Done. End time: $(date)"
echo "  Modified : ${success} / ${total} VCF files"
echo "================================================"
