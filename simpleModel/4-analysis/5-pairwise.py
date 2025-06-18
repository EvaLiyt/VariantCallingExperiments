import os
import re
from Bio import SeqIO
from ete3 import Tree

NUM_SITES = 10000 

# Folders
root_folder = "../1-simulation/root"
true_alignment_folder = "../1-simulation/nexus"
called_alignment_folder = "../2-variantCalling/alignment"
output_folder = "./stats"

def get_sequence_from_nexus(nexus_file, cell_name):
    sequence = None
    in_matrix = False
    with open(nexus_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.lower() == "matrix":
                in_matrix = True
                continue
            if in_matrix:
                if line == ";" or line.lower() == "end;":
                    break
                parts = line.split()
                if len(parts) >= 2:
                    taxon = parts[0]
                    seq = "".join(parts[1:])
                    if taxon == cell_name:
                        seq_clean = seq.replace(",", "").replace(".", "").replace(";", "").replace(" ", "").upper()
                        sequence = seq_clean
                        break
    return sequence

def get_sequence_from_fasta(fasta_file, cell_name):
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id == cell_name:
            return str(record.seq).upper()
    return None

def get_all_names_from_fasta(fasta_file):
    return [record.id for record in SeqIO.parse(fasta_file, "fasta")]

def compute_pairwise_diff(seq1, seq2):
    return sum(1 for i in range(min(len(seq1), len(seq2))) if seq1[i] != seq2[i])

def process_r_folder(r_folder_name):
    results = []

    called_alignment_path = os.path.join(called_alignment_folder, f"alignment_{r_folder_name}.fasta")
    ref_seq_path = os.path.join(root_folder, f"sequenceGeneration_{r_folder_name}_D.nexus")
    true_seq_path = os.path.join(true_alignment_folder, f"sequenceGeneration_{r_folder_name}_A.nexus")

    cell_names = get_all_names_from_fasta(called_alignment_path)

    for cell_name in cell_names:
        called_seq = get_sequence_from_fasta(called_alignment_path, cell_name)
        ref_seq = get_sequence_from_nexus(ref_seq_path, "chr1")  # assuming reference seq named "chr1"
        true_seq = get_sequence_from_nexus(true_seq_path, cell_name)

        if None in [called_seq, ref_seq, true_seq]:
            print(f"[WARNING] Missing sequences for cell {cell_name} in {r_folder_name}")
            continue

        pairwise_diff = compute_pairwise_diff(ref_seq, true_seq)

        results.append((r_folder_name, cell_name, pairwise_diff))

    return results

def main():
    all_results = []
    for i in range(100):
        r_folder_name = f"r{i}"
        print(f"Processing {r_folder_name}...")
        try:
            res = process_r_folder(r_folder_name)
            all_results.extend(res)
        except Exception as e:
            print(f"[ERROR] Processing {r_folder_name}: {e}")

    output_path = os.path.join(output_folder, "pairwise.tsv")
    with open(output_path, 'w') as out_f:
        out_f.write("r_folder\tcell_name\tpairwiseDiff2\n")
        for r, cell, diff in all_results:
            out_f.write(f"{r}\t{cell}\t{diff}\n")

    print(f"Output saved to {output_path}")

if __name__ == "__main__":
    main()
