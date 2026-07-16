#!/usr/bin/env python3

# ============================================================
# generateFasta.py
# Run from: 2-variantCalling/
#
# For each DxVx/merge/variants_r*_dx_vx.vcf:
#   - Load the matching reference sequence from
#     ../1-simulation/root/sequenceGeneration_<r_part>_D.nexus
#   - Convert the reference sequence using unphased rules
#   - Write one FASTA entry per sample column
# Output: DxVx/merge/root_ref_<suffix>.fasta
# ============================================================

import os


def parse_nexus_get_first_seq(nexus_file):
    with open(nexus_file, "r") as f:
        lines = f.readlines()

    in_matrix = False
    sequences = []

    for line in lines:
        line = line.strip()
        if line.lower() == "matrix":
            in_matrix = True
            continue
        if in_matrix:
            if line == ";":
                break
            if line == "":
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            taxon = parts[0]
            seq   = "".join(parts[1:])
            for ch in [",", ".", ";", " ", "\t", "\n", "\r"]:
                seq = seq.replace(ch, "")
            sequences.append((taxon, seq))

    if not sequences:
        raise ValueError(f"No sequences found in matrix of {nexus_file}")

    return sequences[0]


def convert_unphased(sequence):
    mapping = {
        '0': '0', '5': '5', 'a': 'a', 'f': 'f',
        '1': 'M', '4': 'M',
        '2': 'R', '8': 'R',
        '3': 'W', 'c': 'W',
        '6': 'S', '9': 'S',
        '7': 'Y', 'd': 'Y',
        'b': 'K', 'e': 'K',
    }
    converted = []
    for char in sequence.lower():
        if char in mapping:
            converted.append(mapping[char])
        else:
            raise ValueError(f"Invalid character '{char}' in sequence")
    return ''.join(converted)


# ── Paths ─────────────────────────────────────────────────────
base_dir       = os.path.dirname(os.path.abspath(__file__))
ref_seq_folder = os.path.join(base_dir, "..", "..", "1-simulation", "root")

fasta_count = 0
warn_count  = 0

# ── Loop over all DxVx folders ────────────────────────────────
for dxvx_name in sorted(os.listdir(base_dir)):
    dxvx_dir = os.path.join(base_dir, dxvx_name)
    if not (os.path.isdir(dxvx_dir) and dxvx_name.startswith('D')):
        continue

    merge_dir = os.path.join(dxvx_dir, "merge")
    if not os.path.isdir(merge_dir):
        print(f"  [SKIP] {dxvx_name}: no merge/ directory")
        continue

    print(f"\nProcessing: {dxvx_name}")

    for file_name in sorted(os.listdir(merge_dir)):
        if not (file_name.startswith("variants_") and
                file_name.endswith(".vcf")):
            continue

        # e.g. variants_r0_d1_v1.vcf -> suffix = r0_d1_v1
        suffix = file_name.replace("variants_", "").replace(".vcf", "")

        # Extract just the r* part (e.g. "r0" from "r0_d1_v1")
        r_part = suffix.split("_")[0]

        # Locate matching nexus file using only the r* part
        nexus_file = os.path.join(
            ref_seq_folder,
            f"sequenceGeneration_{r_part}_D.nexus"
        )
        if not os.path.exists(nexus_file):
            print(f"  [WARN] Nexus file not found for {r_part}: {nexus_file}")
            warn_count += 1
            continue

        # Parse reference sequence
        try:
            reference_taxon, reference_sequence = \
                parse_nexus_get_first_seq(nexus_file)
        except Exception as e:
            print(f"  [ERROR] Parsing {nexus_file}: {e}")
            warn_count += 1
            continue

        ref_sequence = convert_unphased(reference_sequence)
        print(f"  {suffix}: reference '{reference_taxon}' "
              f"length {len(reference_sequence)}")

        # Read sample names from #CHROM line
        vcf_path = os.path.join(merge_dir, file_name)
        sample_names = None
        with open(vcf_path, "r") as vcf_file:
            for line in vcf_file:
                if line.startswith("#CHROM"):
                    columns      = line.strip().split('\t')
                    sample_names = columns[9:]
                    break

        if sample_names is None:
            print(f"  [WARN] No #CHROM line found in {file_name}")
            warn_count += 1
            continue

        # Write FASTA
        output_fasta_path = os.path.join(
            merge_dir, f"root_ref_{suffix}.fasta"
        )
        with open(output_fasta_path, "w") as output_fasta:
            for sample_name in sample_names:
                output_fasta.write(f">{sample_name}\n")
                output_fasta.write(f"{ref_sequence}\n")

        print(f"    -> root_ref_{suffix}.fasta "
              f"({len(sample_names)} sample(s))")
        fasta_count += 1

print(f"\nFinished. {fasta_count} FASTA file(s) written, "
      f"{warn_count} warning(s).")