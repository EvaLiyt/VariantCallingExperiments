import os

def parse_nexus_get_first_seq(nexus_file):
    """
    Parses a Nexus file to extract the first taxon's sequence from the Matrix block.
    Cleans sequence by removing ',', '.', ';' and whitespace.
    Returns a tuple: (taxon_name, cleaned_sequence)
    """
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
            if line == ";":  # end of matrix block
                break
            if line == "":
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            taxon = parts[0]
            seq = "".join(parts[1:])

            # Clean sequence: remove ',', '.', ';' and whitespace
            for ch in [",", ".", ";", " ", "\t", "\n", "\r"]:
                seq = seq.replace(ch, "")

            sequences.append((taxon, seq))

    if not sequences:
        raise ValueError(f"No sequences found in matrix of {nexus_file}")

    return sequences[0]  # first taxon and cleaned sequence


# Set folders
ref_sequence_folder = "../1-simulation/root"
vcf_folder = "./merge"
output_folder = "./merge"


for file_name in os.listdir(vcf_folder):
    if file_name.startswith("variants_") and file_name.endswith(".vcf"):
        suffix = file_name.replace("variants_", "").replace(".vcf", "")
        nexus_file = os.path.join(ref_sequence_folder, f"sequenceGeneration_{suffix}_D.nexus")

        if not os.path.exists(nexus_file):
            print(f"Warning: Nexus file not found for {suffix}: {nexus_file}")
            continue

        try:
            reference_taxon, reference_sequence = parse_nexus_get_first_seq(nexus_file)
        except Exception as e:
            print(f"Error parsing {nexus_file}: {e}")
            continue

        seq_len = len(reference_sequence)
        unknown_sequence = "?" * seq_len  # <-- this is what you want

        print(f"Reference length for {suffix} ({reference_taxon}): {seq_len}")

        vcf_path = os.path.join(vcf_folder, file_name)
        with open(vcf_path, "r") as vcf_file:
            for line in vcf_file:
                if line.startswith("#CHROM"):
                    columns = line.strip().split('\t')
                    sample_names = columns[9:]  # sample names start from 10th column
                    break
            else:
                print(f"Warning: No #CHROM line found in {vcf_path}")
                continue

        output_fasta_path = os.path.join(output_folder, f"root_ref_{suffix}.fasta")
        with open(output_fasta_path, "w") as output_fasta:
            for sample_name in sample_names:
                output_fasta.write(f">{sample_name}\n")
                output_fasta.write(f"{unknown_sequence}\n")

        print(f"FASTA written: {output_fasta_path}")

print("All FASTA files generated.")
