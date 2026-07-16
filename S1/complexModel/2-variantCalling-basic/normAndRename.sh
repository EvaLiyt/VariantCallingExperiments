#!/bin/bash

# ============================================================
# 7-normAndRename.sh
# Run from: 2-variantCalling/
#
# For each DxVx/modifiedVCF/rx/*.vcf:
#   Step A (bcftools): compress → bgzip, normalise → .norm.vcf.gz, index
#   Step B (Python):   read .norm.vcf.gz, rename #CHROM sample columns
#                      to the base filename, compress+sort+index
#   Output: DxVx/changedVCF/rx/<same base>.norm.vcf.gz  (+.tbi)
# ============================================================

BASE_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "================================================"
echo "Working directory: ${BASE_DIR}"
echo "Start time: $(date)"
echo "================================================"

# ── Load modules ─────────────────────────────────────────────
module load GCCcore/11.3.0
module load BCFtools/1.19-GCC-11.3.0

# ── Step A: compress and normalise ───────────────────────────
echo ""
echo "[A] bcftools: compress and normalise modifiedVCF..."

for dxvx_dir in "${BASE_DIR}"/D*V*/; do
    [[ -d "${dxvx_dir}" ]] || continue
    sample=$(basename "${dxvx_dir}")

    mod_dir="${dxvx_dir}modifiedVCF"
    if [[ ! -d "${mod_dir}" ]]; then
        echo "  [SKIP] ${sample}: no modifiedVCF/ directory"
        continue
    fi

    echo ""
    echo "  ${sample}"

    for rx_dir in "${mod_dir}"/r*/; do
        [[ -d "${rx_dir}" ]] || continue
        rep=$(basename "${rx_dir}")

        for vcf_file in "${rx_dir}"/*.vcf; do
            [[ -f "${vcf_file}" ]] || continue

            fname=$(basename "${vcf_file}")
            base="${fname%.vcf}"

            gz_file="${rx_dir}${base}.vcf.gz"
            norm_file="${rx_dir}${base}.norm.vcf.gz"

            echo " ${rep}/${fname}"

            # Compress
            bcftools view "${vcf_file}" -Oz -o "${gz_file}"

            # Normalise
            bcftools norm -Oz -f ~/data/10K.fna -o "${norm_file}" "${gz_file}"

            # Index
            bcftools index "${norm_file}"

            echo " -> ${base}.norm.vcf.gz (.csi)"
        done
    done
done

echo ""
echo "[A] Normalisation complete."

# ── Step B: Python — rename sample columns → changedVCF ──────
echo ""
echo "[B] Python: renaming sample headers and writing changedVCF..."

python "${BASE_DIR}/renameVCF.py" "${BASE_DIR}"

echo ""
echo "================================================"
echo "Done. End time: $(date)"
echo "================================================"
