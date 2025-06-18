import os

input_folder = "./files"
output_folder = "./comparison"
params = "../1-simulation/pileup"

# Ensure output folder exists
os.makedirs(output_folder, exist_ok=True)

def read_sequences_from_file(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    return [line.strip() for line in lines[1:]]  # skip header line, strip newlines

def compare_sequences(seq1, seq2, seq3):
    differing_sites = []
    length = min(len(seq1), len(seq2), len(seq3))
    for i in range(length):
        if seq1[i] != seq2[i] or seq2[i] != seq3[i] or seq1[i] != seq3[i]:
            differing_sites.append((i+1, seq1[i], seq2[i], seq3[i]))
    return differing_sites

def parse_params(file_name):
    r = os.path.basename(file_name).replace("parameters_", "").replace(".txt", "")
    minD = None
    minVar = None
    try:
        with open(file_name, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith("min_depth:"):
                    minD = line.split(":")[1].strip()
                elif line.startswith("min_var:"):
                    minVar = line.split(":")[1].strip()
    except FileNotFoundError:
        pass
    return r, minD, minVar

def process_r_folder(r_folder_name):
    sub_input_folder = os.path.join(input_folder, r_folder_name)
    sub_output_folder = os.path.join(output_folder, r_folder_name)
    os.makedirs(sub_output_folder, exist_ok=True)
    
    params_file = os.path.join(params, r_folder_name, f"parameters_{r_folder_name}.txt")
    r, minD, minVar = parse_params(params_file)
    
    for fasta_file in os.listdir(sub_input_folder):
        fasta_path = os.path.join(sub_input_folder, fasta_file)
        sequences = read_sequences_from_file(fasta_path)
        if len(sequences) != 3:
            print(f"Skipping '{fasta_file}': requires exactly 3 sequences.")
            continue

        differing_sites = compare_sequences(*sequences)
        cell_name = fasta_file.split(".")[0]
        output_file = os.path.join(sub_output_folder, f"{cell_name}.tsv")
        
        with open(output_file, 'w') as out:
            out.write("site\tref\ttrue\tfinal\n")
            for site in differing_sites:
                out.write(f"{site[0]}\t{site[1]}\t{site[2]}\t{site[3]}\n")

        print(f"Wrote differences for '{fasta_file}' to '{output_file}'.")

def main():
    for i in range(100):
        process_r_folder(f"r{i}")

if __name__ == "__main__":
    main()
