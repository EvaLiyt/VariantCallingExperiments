import os

# Set file paths
reference_sequence_file = os.path.expanduser("~/data/chr1.fna")
vcf_folder = "./merge"
output_folder = "./merge"

# Read the reference FASTA
with open(reference_sequence_file, "r") as ref_file:
    ref_lines = ref_file.readlines()

# Extract reference sequence (excluding header)
reference_sequence = ''.join(line.strip() for line in ref_lines if not line.startswith(">"))

if not reference_sequence:
    print("Warning: Reference sequence is empty.")
    exit(1)
else:
    print(f"Reference sequence loaded successfully with length: {len(reference_sequence)}")

# Ensure output folder exists
os.makedirs(output_folder, exist_ok=True)

# Loop through VCF files like variants_r0.vcf, variants_r1.vcf, etc.
for file_name in os.listdir(vcf_folder):
    if not file_name.startswith("variants_") or not file_name.endswith(".vcf"):
        continue

    file_path = os.path.join(vcf_folder, file_name)
    rx = file_name.split('.')[0].split('_')[1]  # e.g., r0, r1, etc.

    sample_names = []
    with open(file_path, "r") as vcf_file:
        for line in vcf_file:
            if line.startswith("#CHROM"):
                columns = line.strip().split('\t')
                sample_names = columns[9:]  # Columns after FORMAT
                break

    if not sample_names:
        print(f"Warning: No sample names found in {file_name}")
        continue

    output_fasta_path = os.path.join(output_folder, f"new_ref_{rx}.fasta")
    with open(output_fasta_path, "w") as fasta_file:
        for sample in sample_names:
            fasta_file.write(f">{sample}\n")
            fasta_file.write(f"{reference_sequence}\n")

    print(f"Generated FASTA: {output_fasta_path}")

print("All FASTA files generated.")
