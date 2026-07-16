#!/bin/bash
# ============================================================
# 5-filterVCF.sh
# Run from: 2-variantCalling/
#
# 1. Remove *.out and *.log files inside each DxVx/ folder
# 2. For each DxVx/vcResults/rx/*.vcf:
#    - Keep heterozygous calls (0/1) flagged True
#    - Keep SNPs only (ALT length == 1)
#    - Keep FILTER column == '.' (no filter applied)
#    - Keep sites with read depth (AD[ref]+AD[alt]) >= 20
#    Output: same path with .filtered.vcf suffix
# ============================================================
BASE_DIR="$(cd "$(dirname "$0")" && pwd)"
echo "================================================"
echo "Working directory: ${BASE_DIR}"
echo "Start time: $(date)"
echo "================================================"
# ── 1. Clean *.out and *.log files from DxVx folders ────────
echo ""
echo "[1/2] Removing *.out and *.log files..."
for dxvx_dir in "${BASE_DIR}"/D*V*/; do
    [[ -d "${dxvx_dir}" ]] || continue
    sample=$(basename "${dxvx_dir}")
    out_count=$(find "${dxvx_dir}" -maxdepth 1 -name "*.out" | wc -l)
    log_count=$(find "${dxvx_dir}" -maxdepth 1 -name "*.log" | wc -l)
    if [[ $((out_count + log_count)) -gt 0 ]]; then
        find "${dxvx_dir}" -maxdepth 1 -name "*.out" -delete
        find "${dxvx_dir}" -maxdepth 1 -name "*.log" -delete
        echo "  ${sample}: removed ${out_count} *.out  ${log_count} *.log"
    else
        echo "  ${sample}: nothing to clean"
    fi
done
# ── 2. Filter VCF files ──────────────────────────────────────
echo ""
echo "[2/2] Filtering VCF files..."
total=0
success=0
skipped=0
for dxvx_dir in "${BASE_DIR}"/D*V*/; do
    [[ -d "${dxvx_dir}" ]] || continue
    sample=$(basename "${dxvx_dir}")
    vcResults="${dxvx_dir}vcResults"
    if [[ ! -d "${vcResults}" ]]; then
        echo "  [SKIP] ${sample}: no vcResults/ directory"
        continue
    fi
    for rx_dir in "${vcResults}"/r*/; do
        [[ -d "${rx_dir}" ]] || continue
        rep=$(basename "${rx_dir}")
        # Mirror output path: DxVx/filtered/rx/
        out_dir="${dxvx_dir}filtered/${rep}"
        mkdir -p "${out_dir}"
        while IFS= read -r -d '' input; do
            fname=$(basename "${input}")   # keep exact same filename
            o1="${out_dir}/${fname}"
            total=$((total + 1))
            echo "  Processing: ${sample}/vcResults/${rep}/${fname}"
            # Write header lines first, then append filtered variant lines
            grep '^#' "${input}" > "${o1}"
            grep 'True' "${input}" | awk '$7 == "." && length($5) == 1' >> "${o1}"
            line_count=$(wc -l < "${o1}")
            echo "    -> filtered/${rep}/${fname}  (${line_count} lines: headers + retained variants)"
            success=$((success + 1))
        done < <(find "${rx_dir}" -maxdepth 1 -name "*.vcf" -print0)
    done
done
# ── Summary ──────────────────────────────────────────────────
echo ""
echo "================================================"
echo "Done. End time: $(date)"
echo "  VCF files processed : ${success} / ${total}"
echo "================================================"
