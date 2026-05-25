#!/bin/bash
# ============================================================
# 5-filterVCF.sh
# Run from the directory containing step1/
#
# 1. Remove *.out and *.log in current directory
# 2. For each step1/rx/*.sccaller.vcf:
#    - Collect variant sites (GT == 0/1 or 1/1) with SO == True
#    - Output 0-based BED to step2/rx/*_varPos.bed
# ============================================================

echo "Working directory: $(pwd)"
echo "Start time: $(date)"

# ── 1. Clean *.out and *.log ─────────────────────────────────
echo ""
echo "[1/2] Removing *.out and *.log files..."
out_count=$(find . -maxdepth 1 -name "*.out" | wc -l)
log_count=$(find . -maxdepth 1 -name "*.log" | wc -l)
find . -maxdepth 1 -name "*.out" -delete
find . -maxdepth 1 -name "*.log" -delete
echo "  Removed ${out_count} *.out and ${log_count} *.log"

# ── 2. Generate per-cell BED files ───────────────────────────
echo ""
echo "[2/2] Generating BED files..."
total=0
success=0

for rx_dir in step1/r*/; do
    [[ -d "$rx_dir" ]] || continue
    rep=$(basename "$rx_dir")
    out_dir="step2/${rep}"
    mkdir -p "$out_dir"

    for vcf in "${rx_dir}"*.sccaller.vcf; do
        [[ -f "$vcf" ]] || continue
        fname=$(basename "$vcf" .sccaller.vcf)
        out_bed="${out_dir}/${fname}_varPos.bed"
        total=$((total + 1))

        echo "  Processing: $vcf"
        grep -v '^#' "$vcf" | awk '
        {
            split($10, s, ":")
            if ((s[1] == "0/1" || s[1] == "1/1") && s[2] == "True")
                print $1 "\t" ($2 - 1) "\t" $2
        }' | sort -k1,1 -k2,2n | uniq > "$out_bed"

        lines=$(wc -l < "$out_bed")
        echo "    -> ${out_bed} (${lines} variant sites)"
        success=$((success + 1))
    done
done

# ── Summary ──────────────────────────────────────────────────
echo ""
echo "================================================"
echo "Done. End time: $(date)"
echo "  BED files written: ${success} / ${total}"
echo "================================================"