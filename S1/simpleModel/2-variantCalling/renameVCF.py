#!/usr/bin/env python3

# ============================================================
# renameVCF.py
# Run from: 2-variantCalling/
# Called by: 7-normAndRename.sh
#
# For each DxVx/modifiedVCF/rx/*.norm.vcf.gz:
#   - Rename #CHROM sample columns to the base filename
#   - bgzip compress, sort, tabix index
# Output: DxVx/changedVCF/rx/<same filename> + .tbi
# ============================================================

import os, sys, gzip, shutil, subprocess

base_dir = sys.argv[1]

for dxvx_name in sorted(os.listdir(base_dir)):
    dxvx_dir = os.path.join(base_dir, dxvx_name)
    if not (os.path.isdir(dxvx_dir) and dxvx_name.startswith('D')):
        continue

    mod_dir = os.path.join(dxvx_dir, "modifiedVCF")
    if not os.path.isdir(mod_dir):
        print(f"  [SKIP] {dxvx_name}: no modifiedVCF/")
        continue

    changed_base = os.path.join(dxvx_dir, "changedVCF")
    if os.path.exists(changed_base):
        shutil.rmtree(changed_base)
    os.makedirs(changed_base)
    print(f"\n  {dxvx_name}")

    for rep in sorted(os.listdir(mod_dir)):
        rep_dir = os.path.join(mod_dir, rep)
        if not os.path.isdir(rep_dir):
            continue

        out_rep_dir = os.path.join(changed_base, rep)
        os.makedirs(out_rep_dir, exist_ok=True)

        for filename in sorted(os.listdir(rep_dir)):
            if not filename.endswith(".norm.vcf.gz"):
                continue

            input_path = os.path.join(rep_dir, filename)
            base_name  = filename.replace(".norm.vcf.gz", "")
            out_vcf    = os.path.join(out_rep_dir, filename.replace(".norm.vcf.gz", ".norm.vcf"))  # FIXED
            out_gz     = os.path.join(out_rep_dir, filename)   # .norm.vcf.gz
            out_tbi    = out_gz + ".tbi"

            print(f"    {rep}/{filename}")

            # Rename #CHROM sample columns to base_name
            with gzip.open(input_path, 'rt') as infile, \
                 open(out_vcf, 'w') as outfile:
                for line in infile:
                    if line.startswith("#CHROM"):
                        headers = line.strip().split('\t')
                        headers[9:] = [base_name] * len(headers[9:])
                        line = '\t'.join(headers) + '\n'
                    outfile.write(line)

            # bgzip compress
            subprocess.run(f"bgzip -f {out_vcf}", shell=True, check=True,
                           executable='/bin/bash')
            # sort in-place
            tmp_sorted = out_gz + ".tmp.gz"
            subprocess.run(f"bcftools sort -Oz -o {tmp_sorted} {out_gz}",
                           shell=True, check=True, executable='/bin/bash')
            os.replace(tmp_sorted, out_gz)
            # tabix index
            subprocess.run(f"tabix -p vcf {out_gz}",
                           shell=True, check=True, executable='/bin/bash')

            print(f"      -> changedVCF/{rep}/{filename} (.tbi)")

print("\nFinished renaming, compressing, and indexing VCF files.")
