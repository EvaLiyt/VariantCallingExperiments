import os
import re
from Bio import SeqIO

NUM_SITES = 10000 

# Folders
root_folder = "../../1-simulation/root"
true_alignment_folder = "../../1-simulation/nexus"
called_alignment_folder = "../../2-variantCalling/filterMethods/method2/alignment"
output_folder = "../summary"

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

def get_all_sequences_from_nexus(nexus_file):
    sequences = {}
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
                    seq_clean = seq.replace(",", "").replace(".", "").replace(";", "").replace(" ", "").upper()
                    sequences[taxon] = seq_clean
    return sequences

def compute_hammingDist_diff(seq1, seq2):
    return sum(1 for i in range(min(len(seq1), len(seq2))) if seq1[i] != seq2[i])

def process_r_folder(r_folder_name):
    results = []

    true_seq_path = os.path.join(true_alignment_folder, f"sequenceGeneration_{r_folder_name}_A.nexus")

    # Load all sequences from nexus
    sequences = get_all_sequences_from_nexus(true_seq_path)
    cell_names = list(sequences.keys())

    # Pairwise comparisons
    for i in range(len(cell_names)):
        for j in range(i + 1, len(cell_names)):
            cell1, cell2 = cell_names[i], cell_names[j]
            seq1, seq2 = sequences[cell1], sequences[cell2]
            hammingDist_diff = compute_hammingDist_diff(seq1, seq2)
            results.append((r_folder_name, cell1, cell2, hammingDist_diff))

    return results

def main():
    os.makedirs(output_folder, exist_ok=True)
    all_results = []
    for i in range(100):
        r_folder_name = f"r{i}"
        print(f"Processing {r_folder_name}...")
        try:
            res = process_r_folder(r_folder_name)
            all_results.extend(res)
        except Exception as e:
            print(f"[ERROR] Processing {r_folder_name}: {e}")

    output_path = os.path.join(output_folder, "hammingDist_method2.tsv")
    with open(output_path, 'w') as out_f:
        out_f.write("r_folder\tcell1\tcell2\thammingDist\n")
        for r, c1, c2, diff in all_results:
            out_f.write(f"{r}\t{c1}\t{c2}\t{diff}\n")

    print(f"Output saved to {output_path}")

if __name__ == "__main__":
    main()
