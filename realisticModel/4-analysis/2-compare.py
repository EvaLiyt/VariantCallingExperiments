import os
import shutil

input_folder = "./files"
output_folder = "./comparison"
site_folder = "./sites"
params = "../1-simulation/pileup"
vcf_folder = "../2-variantCalling/filtered"

print("Current working directory:", os.getcwd())

# Ensure output directories exist and clean them
for folder in [output_folder, site_folder]:
    if not os.path.exists(folder):
        os.makedirs(folder)
        print(f"Directory '{folder}' created.")
    else:
        # Remove all files and subfolders
        for filename in os.listdir(folder):
            file_path = os.path.join(folder, filename)
            try:
                if os.path.isfile(file_path):
                    os.remove(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print(f"Error removing {file_path}: {e}")


def read_sequences_from_file(file_path):
    """Reads sequences from a file and returns them as a list."""
    with open(file_path, 'r') as file:
        sequences = file.readlines()
        # Skip the first line (assuming it's a header) and strip newlines
        return [seq.strip() for seq in sequences[1:]]  # Assuming sequences start from the second line


def compare_sequences(seq1, seq2, seq3):
    """Compare true and final sequences and find differing sites."""
    differing_sites = []
    for i in range(min(len(seq1), len(seq2), len(seq3))):  # Compare only up to the shortest sequence
        if seq1[i] != seq2[i] or seq2[i] != seq3[i]:
            differing_sites.append((i + 1, seq1[i], seq2[i], seq3[i]))  # Store 1-based index and the differing bases
    return differing_sites


def parse_params(file_name):
    base_name = os.path.basename(file_name)
    if not base_name.startswith("parameters_") or not base_name.endswith(".txt"):
        print(f"File name '{file_name}' does not match expected format.")
        return None, None, None
    r = base_name.replace("parameters_", "").replace(".txt", "")
    minD = None
    minVar = None

    try:
        with open(file_name, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith("min_depth:"):
                    minD = line.split(":")[1].strip()
                elif line.startswith("min_var:"):
                    minVar = line.split(":")[1].strip()
    except FileNotFoundError:
        print(f"File '{file_name}' not found.")
        return r, None, None

    print(f"Parsed r: {r}, minD: {minD}, minVar: {minVar}")
    return r, minD, minVar


def process_r_folder(r_folder_name):
    """Process the r folder and generate outputs in respective subfolders."""
    sub_input_folder = os.path.join(input_folder, r_folder_name)  # Access input files in the subfolder
    sub_output_folder = os.path.join(output_folder, r_folder_name)
    sub_site_folder = os.path.join(site_folder, r_folder_name)

    # Create subdirectories if they don't exist
    os.makedirs(sub_output_folder, exist_ok=True)
    os.makedirs(sub_site_folder, exist_ok=True)

    params_file = os.path.join(params, r_folder_name, f"parameters_{r_folder_name}.txt")
    r, minD, minVar = parse_params(params_file)
    for fasta_file in os.listdir(sub_input_folder):
        fasta_path = os.path.join(sub_input_folder, fasta_file)

        sequences = read_sequences_from_file(fasta_path)
        if len(sequences) != 3:
            print(f"File '{fasta_file}' does not contain exactly 3 sequences. Skipping.")
            continue

        # Compare the three sequences
        differing_sites = compare_sequences(sequences[0], sequences[1], sequences[2])

        # Write differing sites to output file
        cellName = fasta_file.split(".")[0]
        output_file_path = os.path.join(sub_output_folder, f"{cellName}.tsv")
        with open(output_file_path, 'w') as output_file:
            output_file.write("site\tref\ttrue\tfinal\n")
            for site in differing_sites:
                output_file.write(f"{site[0]}\t{site[1]}\t{site[2]}\t{site[3]}\n")

        print(f"Differing sites for '{fasta_file}' written to '{output_file_path}'.")

    # Process VCF files in the VCF folder
    vcf_r_folder = os.path.join(vcf_folder, r_folder_name)
    if not os.path.exists(vcf_r_folder):
        print(f"VCF folder for '{r_folder_name}' not found. Skipping.")
        return

    for vcf_file in os.listdir(vcf_r_folder):
        vcf_path = os.path.join(vcf_r_folder, vcf_file)

        pos_list_positive = []
        pos_list_negative = []
        with open(vcf_path, 'r') as vcf:
            for line in vcf:
                if line.startswith("#"):
                    continue
                columns = line.strip().split("\t")
                if columns[9].split(":")[0] == "0/0":
                    pos_list_negative.append(columns[1])
                else:
                    pos_list_positive.append(columns[1])

        # Write POS data to site folder
        parts = vcf_file.split(".")[0].split("_")
        name = f"{'_'.join(parts[1:])}_{parts[0]}_minD{minD}_minVar{minVar}"

        site_file_path = os.path.join(sub_site_folder, f"{name}.tsv")
        with open(site_file_path, 'w') as site_file:
            site_file.write("positive_sites\n")
            site_file.write("\n".join(pos_list_positive) + "\n")
            site_file.write("negative_sites\n")
            site_file.write("\n".join(pos_list_negative) + "\n")

        print(f"VCF sites for '{vcf_file}' written to '{site_file_path}'.")


def main():
    # Process r0 to r99
    for i in range(100):
        r_folder_name = f"r{i}"
        process_r_folder(r_folder_name)


if __name__ == "__main__":
    main()
