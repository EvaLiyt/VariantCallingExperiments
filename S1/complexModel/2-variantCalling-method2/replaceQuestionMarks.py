#!/usr/bin/env python3
# ============================================================
# alignFasta.py
# Run from: 2-variantCalling/
#
# For each DxVx/merge/root_ref_<suffix>.fasta:
#   - Load reference sequences
#   - Replace '?' in raw_alignment_<suffix>.fasta with ref base
#   - Rename sequence IDs
# Output: DxVx/alignment/alignment_<suffix>.fasta
# ============================================================
import os
import shutil
from Bio import SeqIO
from Bio.Seq import Seq

base_dir = os.path.dirname(os.path.abspath(__file__))
saved_count = 0
warn_count  = 0

for dxvx_name in sorted(os.listdir(base_dir)):
    dxvx_dir = os.path.join(base_dir, dxvx_name)
    if not (os.path.isdir(dxvx_dir) and dxvx_name.startswith('D')):
        continue
    merge_dir  = os.path.join(dxvx_dir, "merge")
    output_dir = os.path.join(dxvx_dir, "alignment")
    if not os.path.isdir(merge_dir):
        print(f"  [SKIP] {dxvx_name}: no merge/ directory")
        continue
    print(f"\nProcessing: {dxvx_name}")
    # Reset alignment folder for this sample
    if os.path.exists(output_dir):
        print(f"  Deleting existing {dxvx_name}/alignment/")
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)
    print(f"  Created {dxvx_name}/alignment/")
    ref_files = [f for f in os.listdir(merge_dir)
                 if f.startswith("root_ref") and f.endswith(".fasta")]
    for ref_fasta in sorted(ref_files):
        # root_ref_r0_d1_v1.fasta -> suffix = r0_d1_v1
        suffix = '_'.join(ref_fasta.split('_')[2:]).split('.')[0]
        ref_fasta_path   = os.path.join(merge_dir, ref_fasta)
        input_fasta_path = os.path.join(merge_dir,
                                        f"raw_alignment_{suffix}.fasta")
        output_fasta_path = os.path.join(output_dir,
                                         f"alignment_{suffix}.fasta")
        if not os.path.exists(input_fasta_path):
            print(f"  [WARN] raw_alignment_{suffix}.fasta not found")
            warn_count += 1
            continue
        print(f"  {suffix}: raw_alignment -> alignment_{suffix}.fasta")
        # Load reference sequences into memory
        ref_dict = {rec.id: str(rec.seq)
                    for rec in SeqIO.parse(ref_fasta_path, "fasta")}
        # Process alignment file
        modified_records = []
        for rec in SeqIO.parse(input_fasta_path, "fasta"):
            ref_seq = ref_dict.get(rec.id)
            if ref_seq:
                # Replace '?' with the corresponding reference base
                seq_list = list(str(rec.seq))
                for i, base in enumerate(seq_list):
                    if base == '?' and i < len(ref_seq):
                        seq_list[i] = ref_seq[i]
                new_seq = ''.join(seq_list)
            else:
                new_seq = str(rec.seq)   # no replacement if ref missing
            # Rename ID
            if 'clade' in rec.id:
                rec.id = '_'.join(rec.id.split('_')[:2])
            else:
                rec.id = rec.id.split('_')[0]
            rec.description = ""
            rec.seq = Seq(new_seq)
            modified_records.append(rec)
        # Write output
        with open(output_fasta_path, "w") as out:
            for rec in modified_records:
                out.write(f">{rec.id}\n{rec.seq}\n")
        print(f"    Saved: alignment_{suffix}.fasta")
        saved_count += 1
print(f"\nFinished. {saved_count} alignment file(s) saved, "
      f"{warn_count} warning(s).")