import os
import re
import csv
from collections import defaultdict

# Directories
vcf_directory = "../2-variantCalling/filtered"
param_directory = "../1-simulation/pileup"
output_directory = "./stats"
root_dir = "../1-simulation/root"    
nexus_dir = "../1-simulation/nexus" 
output_csv = os.path.join(output_directory, "var_counts.csv")

print("Starting to gather results...\n")
def read_nexus_plaintext(filepath):
    """
    Manually parse a NEXUS file to extract sequences.
    Assumes sequences are in MATRIX block, one per line:
    seq_id  sequence
    Ends reading at ';'
    """
    sequences = {}
    with open(filepath) as f:
        in_matrix = False
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.lower().startswith("begin") or line.lower().startswith("end"):
                continue
            if line.lower() == "matrix":
                in_matrix = True
                continue
            if in_matrix:
                if line == ";":
                    break
                parts = line.split()
                if len(parts) < 2:
                    continue
                seq_id = parts[0]
                seq = "".join(parts[1:])
                # Remove unwanted characters , . ;
                seq = seq.replace(",", "").replace(".", "").replace(";", "")
                sequences[seq_id] = seq
    return sequences

def hamming_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of the same length")
    return sum(el1 != el2 for el1, el2 in zip(seq1, seq2))

def compare_sequences(reference_str, sequences):
    ref_length = len(reference_str)
    differences = {}
    total_different_sites = 0
    for seq_id, sequence in sequences.items():
        if len(sequence) != ref_length:
            raise ValueError(
                f"Sequence length mismatch: {seq_id} length {len(sequence)}, reference length {ref_length}")
        diff_count = hamming_distance(reference_str, sequence)
        differences[seq_id] = diff_count
        total_different_sites += diff_count
    return differences, total_different_sites

def process_r_value(r, root_dir="../1-simulation/root", nexus_dir="../1-simulation/nexus"):
    root_file = os.path.join(root_dir, f"sequenceGeneration_r{r}_D.nexus")
    nexus_file = os.path.join(nexus_dir, f"sequenceGeneration_r{r}_A.nexus")

    if not os.path.exists(root_file) or not os.path.exists(nexus_file):
        print(f"Files for r={r} missing, skipping.")
        return None

    root_sequences = read_nexus_plaintext(root_file)
    nexus_sequences = read_nexus_plaintext(nexus_file)

    reference_id, reference_seq = next(iter(root_sequences.items()))

    differences, total_diff_sites = compare_sequences(reference_seq, nexus_sequences)

    # Separate differences by clade prefix
    clade_diffs = [diff for seq_id, diff in differences.items() if seq_id.startswith("clade")]
    non_clade_diffs = [diff for seq_id, diff in differences.items() if not seq_id.startswith("clade")]

    clade_mean = np.mean(clade_diffs) if clade_diffs else 0
    clade_std = np.std(clade_diffs) if clade_diffs else 0
    non_clade_mean = np.mean(non_clade_diffs) if non_clade_diffs else 0
    non_clade_std = np.std(non_clade_diffs) if non_clade_diffs else 0
    avg_diff_sites = total_diff_sites / len(nexus_sequences) if nexus_sequences else 0

    return {
        "r_value": r,
        "clade_mean": clade_mean,
        "clade_std": clade_std,
        "non_clade_mean": non_clade_mean,
        "non_clade_std": non_clade_std,
        "avg_diff_sites": avg_diff_sites
    }

def read_parameters(r_value):
    param_file_path = os.path.join(param_directory, f"r{r_value}", f"parameters_r{r_value}.txt")
    if not os.path.exists(param_file_path):
        print(f"Parameter file not found for r{r_value}: {param_file_path}")
        return None, None
    with open(param_file_path) as f:
        lines = f.readlines()
        minD = int(lines[0].split(":")[1].strip())
        minVar = int(lines[1].split(":")[1].strip())
    return str(minD), str(minVar)

def clean_sequence(seq):
    return ''.join(ch for ch in seq if ch not in {',', '.', ';'})

def count_mutations(ref, query):
    clean_ref = clean_sequence(ref)
    clean_query = clean_sequence(query)
    length = min(len(clean_ref), len(clean_query))
    return sum(1 for i in range(length) if clean_ref[i] != clean_query[i])


def read_custom_fasta(filepath):
    with open(filepath, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
        if not lines or not lines[0].startswith('>'):
            raise ValueError(f"File {filepath} does not start with '>' header")
        if len(lines) < 4:
            raise ValueError(f"File {filepath} has less than 3 sequence lines")
        return lines[0][1:], lines[1:4]

# Step 1: mutation stats
mutation_stats = defaultdict(lambda: {'normal': {}, 'tumour': {}})

# Step 2: variant stats (store values per file to average later)
variant_stats = defaultdict(lambda: {
    "normal_total_variants": [],
    "normal_zero_zero": [],
    "normal_non_zero": [],
    "tumour_total_variants": [],
    "tumour_zero_zero": [],
    "tumour_non_zero": []
})

for root, _, files in os.walk(vcf_directory):
    print(f"Processing folder: {root}")
    for filename in files:
        if not filename.endswith(".vcf"):
            continue

        match = re.match(r"r(\d+)_((clade_)?.+)\.vcf", filename)
        if not match:
            print(f"Skipping unexpected filename: {filename}")
            continue

        r_value = match.group(1)
        sample_name = match.group(2)
        is_tumour = sample_name.startswith("clade_")

        minD_value, minVar_value = read_parameters(r_value)
        if minD_value is None or minVar_value is None:
            print(f"Skipping r{r_value} due to missing parameters.")
            continue

        filepath = os.path.join(root, filename)
        total = zero_zero = non_zero = 0

        with open(filepath, 'r') as vcf_file:
            for line in vcf_file:
                if line.startswith("#"):
                    continue
                total += 1
                genotype = line.strip().split('\t')[9].split(':')[0]
                if genotype == "0/0":
                    zero_zero += 1
                else:
                    non_zero += 1

        key = (r_value, minD_value, minVar_value)
        stats = variant_stats[key]
        if is_tumour:
            stats["tumour_total_variants"].append(total)
            stats["tumour_zero_zero"].append(zero_zero)
            stats["tumour_non_zero"].append(non_zero)
        else:
            stats["normal_total_variants"].append(total)
            stats["normal_zero_zero"].append(zero_zero)
            stats["normal_non_zero"].append(non_zero)

        print(f"{filename}: r{r_value}, tumour={is_tumour}, total={total}")

# Step 3: Write averaged results
print(f"\nWriting combined results to: {output_csv}")
os.makedirs(os.path.dirname(output_csv), exist_ok=True)

with open(output_csv, 'w', newline='') as out_file:
    fieldnames = [
        "r_value", "minD_value", "minVar_value",
        "avg_total_variants_normal", "avg_zero_zero_normal", "avg_non_zero_normal",
        "avg_total_variants_tumour", "avg_zero_zero_tumour", "avg_non_zero_tumour",
        "avg_mutation_normal", "avg_mutation_tumour"
    ]
    writer = csv.DictWriter(out_file, fieldnames=fieldnames)
    writer.writeheader()

    for (r_value, minD_value, minVar_value), stats in sorted(variant_stats.items(), key=lambda x: (int(x[0][0]), int(x[0][1]), int(x[0][2]))):
        def avg(lst):
            return round(sum(lst) / len(lst), 2) if lst else 0

        mutation_normal = mutation_stats[r_value]['normal'].values()
        mutation_tumour = mutation_stats[r_value]['tumour'].values()

        writer.writerow({
            "r_value": r_value,
            "minD_value": minD_value,
            "minVar_value": minVar_value,
            "avg_total_variants_normal": avg(stats["normal_total_variants"]),
            "avg_zero_zero_normal": avg(stats["normal_zero_zero"]),
            "avg_non_zero_normal": avg(stats["normal_non_zero"]),
            "avg_total_variants_tumour": avg(stats["tumour_total_variants"]),
            "avg_zero_zero_tumour": avg(stats["tumour_zero_zero"]),
            "avg_non_zero_tumour": avg(stats["tumour_non_zero"]),
            "avg_mutation_normal": avg(list(mutation_normal)),
            "avg_mutation_tumour": avg(list(mutation_tumour))
        })
        

