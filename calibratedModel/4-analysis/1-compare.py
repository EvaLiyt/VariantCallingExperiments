import os
import shutil
from Bio import SeqIO

root_folder="../1-simulation/root"
true_alignment_folder="../1-simulation/nexus"
called_alignment_folder="../2-variantCalling/alignment"
params = "../1-simulation/pileup"
vcf_folder = "../2-variantCalling/filtered"

comparison_folder = "./comparison"
site_folder = "./sites"


print("Current working directory:", os.getcwd())

# Ensure output directories exist and clean them
for folder in [comparison_folder, site_folder]:
    if not os.path.exists(folder):
        os.makedirs(folder)
        print(f"Directory '{folder}' created.")
    else:
        for filename in os.listdir(folder):
            file_path = os.path.join(folder, filename)
            try:
                if os.path.isfile(file_path):
                    os.remove(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print(f"Error removing {file_path}: {e}")
                
# get information from files                
def get_sequence_from_nexus(nexus_file, cell_name):
    """
    Reads a Nexus file and returns the sequence for the given cell_name.
    Returns None if cell_name not found.
    """
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
                    # End of matrix block
                    break
                parts = line.split()
                if len(parts) >= 2:
                    taxon = parts[0]
                    seq = "".join(parts[1:])  # In case sequence is split
                    if taxon == cell_name:
                        # Clean sequence if needed, e.g. remove , . ;
                        seq_clean = seq.replace(",", "").replace(".", "").replace(";", "").replace(" ", "").upper()
                        sequence = seq_clean
                        break
    return sequence

def get_sequence_from_fasta(fasta_file, cell_name):
    """
    Reads a FASTA file and returns the sequence string for the given cell_name (sequence ID).
    Returns None if cell_name not found.
    """
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id == cell_name:
            return str(record.seq).upper()
    return None
    
def get_all_names_from_fasta(fasta_file):
    """
    Reads a FASTA file and returns a list of all sequence IDs (cell names).
    """
    names = [record.id for record in SeqIO.parse(fasta_file, "fasta")]
    return names
    

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
  
  
# comparing functions  
def compare_sequences(seq1, seq2, seq3):
    """Compare three sequences and find differing sites."""
    differing_sites = []
    for i in range(min(len(seq1), len(seq2), len(seq3))):
        if seq1[i] != seq2[i] or seq2[i] != seq3[i] or seq1[i] != seq3[i]:
            differing_sites.append((i + 1, seq1[i], seq2[i], seq3[i]))  # 1-based index
    return differing_sites

  

def process_r_folder(r_folder_name):
    """Process the r folder and generate outputs in respective subfolders."""
    sub_comparison_folder = os.path.join(comparison_folder, r_folder_name)
    sub_site_folder = os.path.join(site_folder, r_folder_name)

    os.makedirs(sub_comparison_folder, exist_ok=True)
    os.makedirs(sub_site_folder, exist_ok=True)
    
    params_file = os.path.join(params, r_folder_name, f"parameters_{r_folder_name}.txt")
    r, minD, minVar = parse_params(params_file)

    seq1_seq2_diff_summary = []  # To collect counts for all cells in this r folder
    called_alignment_path = os.path.join(called_alignment_folder, f"alignment_{r_folder_name}.fasta")
    ref_seq_path =  os.path.join(root_folder, f"sequenceGeneration_{r_folder_name}_D.nexus")
    true_seq_path = os.path.join(true_alignment_folder, f"sequenceGeneration_{r_folder_name}_A.nexus")
    
    names = get_all_names_from_fasta(called_alignment_path)
    
    for cell_name in names:
        called_seq = get_sequence_from_fasta(called_alignment_path, cell_name)
        ref_seq = get_sequence_from_nexus(ref_seq_path, "chr1")
        true_seq = get_sequence_from_nexus(true_seq_path, cell_name)

        differing_sites = compare_sequences(ref_seq, true_seq, called_seq)

        output_file_path = os.path.join(sub_comparison_folder, f"{cell_name}.tsv")
        with open(output_file_path, 'w') as output_file:
            output_file.write("site\tref\ttrue\tfinal\n")
            for site in differing_sites:
                output_file.write(f"{site[0]}\t{site[1]}\t{site[2]}\t{site[3]}\n")

        print(f"Differing sites for '{cell_name}_{r_folder_name}' written to '{output_file_path}'.")

        # Count differences between seq1 and seq2
        seq1_seq2_diff_count = sum(
            1 for i in range(min(len(ref_seq), len(true_seq)))
            if ref_seq[i] != true_seq[i]
        )
        seq1_seq2_diff_summary.append((cell_name, seq1_seq2_diff_count))

    # Write summary file for all cells in this r folder
    r_summary_path = os.path.join(sub_comparison_folder, f"{r_folder_name}_pairwise.tsv")
    with open(r_summary_path, 'w') as summary_file:
        summary_file.write("cell_name\tpairwiseDiff2\n")
        for cell_name, count in seq1_seq2_diff_summary:
            summary_file.write(f"{cell_name}\t{count}\n")

    print(f"Summary file for '{r_folder_name}' saved to '{r_summary_path}'.")

    # Process VCF files
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
    for i in range(100):  # Process r0 to r99
        r_folder_name = f"r{i}"
        process_r_folder(r_folder_name)

if __name__ == "__main__":
    main()
