#!/usr/bin/env python3

# ============================================================
# processVariants.py
# Run from: 2-variantCalling/
#
# For each DxVx/merge/merged_r*_dx_vx.vcf:
#   - Revert ##FORMAT=<ID=PL Number=1 back to Number=G
#   - Strip the DxVx prefix from sample column names
#     (e.g. "D1V1_sample" -> "sample")
# Output: DxVx/merge/variants_r*_dx_vx.vcf
# ============================================================

import os

base_dir = os.path.dirname(os.path.abspath(__file__))

processed = 0

for dxvx_name in sorted(os.listdir(base_dir)):
    dxvx_dir = os.path.join(base_dir, dxvx_name)
    if not (os.path.isdir(dxvx_dir) and dxvx_name.startswith('D')):
        continue

    merge_dir = os.path.join(dxvx_dir, "merge")
    if not os.path.isdir(merge_dir):
        print(f"  [SKIP] {dxvx_name}: no merge/ directory")
        continue

    print(f"\nProcessing: {dxvx_name}")

    for filename in sorted(os.listdir(merge_dir)):
        if not (filename.startswith("merged_") and filename.endswith(".vcf")):
            continue

        input_path = os.path.join(merge_dir, filename)

        # merged_r0_d1_v1.vcf -> variants_r0_d1_v1.vcf
        output_filename = filename.replace("merged_", "variants_", 1)
        output_path = os.path.join(merge_dir, output_filename)

        if os.path.exists(output_path):
            print(f"  {output_filename} exists — overwriting...")
            os.remove(output_path)

        print(f"  {filename} -> {output_filename}")

        with open(input_path, 'r') as infile, \
             open(output_path, 'w') as outfile:
            for line in infile:
                # Revert Number=1 back to Number=G in PL FORMAT header
                if line.startswith("##FORMAT=<ID=PL"):
                    line = line.replace("Number=1", "Number=G")

                # Strip DxVx prefix from sample column names
                if line.startswith("#CHROM"):
                    parts  = line.strip().split('\t')
                    header = parts[:9]
                    cells  = [cell.split('_', 1)[1] if '_' in cell else cell
                               for cell in parts[9:]]
                    line = '\t'.join(header + cells) + '\n'

                outfile.write(line)

        processed += 1
        print(f"    Saved: {output_filename}")

print(f"\nFinished. {processed} file(s) processed.")
