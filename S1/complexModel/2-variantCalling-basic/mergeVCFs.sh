#!/bin/bash

# ============================================================
# 8-mergeVCF.sh
# Run from: 2-variantCalling/
#
# For each DxVx/changedVCF/rx/:
#   - Merge all *.norm.vcf.gz files within that replicate folder
# Output: DxVx/merge/merged_r0_d1_v1.vcf  (one per rx per DxVx)
# ============================================================

BASE_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "================================================"
echo "Working directory: ${BASE_DIR}"
echo "Start time: $(date)"
echo "================================================"

# ── Load modules ─────────────────────────────────────────────
module load GCCcore/11.3.0
module load BCFtools/1.19-GCC-11.3.0

merge_count=0
error_count=0

for dxvx_dir in "${BASE_DIR}"/D*V*/; do
    [[ -d "${dxvx_dir}" ]] || continue
    sample=$(basename "${dxvx_dir}")          # e.g. D1V1

    changed_dir="${dxvx_dir}changedVCF"
    if [[ ! -d "${changed_dir}" ]]; then
        echo "  [SKIP] ${sample}: no changedVCF/ directory"
        continue
    fi

    # Parse D1V1 -> d1_v1 for the output filename
    dx_vx=$(echo "${sample}" | sed 's/\(D[0-9]*\)\(V[0-9]*\)/\L\1_\2/')

    # Reset merge folder for this sample
    merge_dir="${dxvx_dir}merge"
    rm -rf "${merge_dir}"
    mkdir -p "${merge_dir}"

    echo ""
    echo "Processing: ${sample}"

    for rx_dir in "${changed_dir}"/r*/; do
        [[ -d "${rx_dir}" ]] || continue
        rep=$(basename "${rx_dir}")           # e.g. r0, r1

        # Collect all .norm.vcf.gz files for this replicate
        vcf_list="${merge_dir}/vcf_list_${rep}.txt"
        > "${vcf_list}"                       # reset list

        for vcf_file in "${rx_dir}"/*.norm.vcf.gz; do
            [[ -f "${vcf_file}" ]] || continue
            echo "${vcf_file}" >> "${vcf_list}"
        done

        file_count=$(wc -l < "${vcf_list}")
        if [[ "${file_count}" -eq 0 ]]; then
            echo "  [SKIP] ${rep}: no .norm.vcf.gz files found"
            rm -f "${vcf_list}"
            continue
        fi

        out_vcf="${merge_dir}/merged_${rep}_${dx_vx}.vcf"
        echo "  ${rep}: merging ${file_count} file(s) -> merged_${rep}_${dx_vx}.vcf"

        bcftools merge --file-list "${vcf_list}" -o "${out_vcf}"

        if [[ $? -eq 0 ]]; then
            echo "    OK"
            merge_count=$((merge_count + 1))
        else
            echo "    ERROR merging ${sample}/${rep}" >&2
            error_count=$((error_count + 1))
        fi

        rm -f "${vcf_list}"
    done
done

echo ""
echo "================================================"
echo "Done. End time: $(date)"
echo "  Merged : ${merge_count}"
echo "  Errors : ${error_count}"
echo "================================================"
