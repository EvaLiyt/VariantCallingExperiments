module load Python

shopt -s nullglob

for folder in ../2-variantCalling-basic-method2/D*V*/; do
    fname=$(basename "$folder")
    [[ "$fname" == "D1V1" ]] && continue

    src_merge="${folder}merge"
    dst_merge="./${fname}/merge"

    if [[ ! -d "${src_merge}" ]]; then
        echo "[WARN] ${fname}: no merge/ in source, skipping" >&2
        continue
    fi

    vcf_files=("${src_merge}"/variants_*.vcf)
    if [[ ${#vcf_files[@]} -eq 0 ]]; then
        echo "[WARN] ${fname}: no variants_*.vcf found in ${src_merge}, skipping" >&2
        continue
    fi

    fasta_files=("${src_merge}"/root_ref_*.fasta)
    if [[ ${#fasta_files[@]} -eq 0 ]]; then
        echo "[WARN] ${fname}: no root_ref_*.fasta found in ${src_merge}, skipping" >&2
        continue
    fi

    mkdir -p "${dst_merge}"
    cp "${vcf_files[@]}" "${dst_merge}/"
    cp "${fasta_files[@]}" "${dst_merge}/"
    echo "[OK] ${fname}: copied ${#vcf_files[@]} VCF file(s)"
done

shopt -u nullglob

chmod +x ./getMappedFile.sh
./getMappedFile.sh