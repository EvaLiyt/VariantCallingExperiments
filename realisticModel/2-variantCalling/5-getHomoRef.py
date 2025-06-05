import os
import gzip

# Set file paths
reference_sequence_file = os.path.expanduser("~/data/chr1_homo.fna")
vcf_folder = "./merge"
output_folder = "./merge"

# Read the reference fasta
with open(reference_sequence_file, "r") as ref_file:
    ref_lines = ref_file.readlines()

# Get the reference sequence (excluding header line)
reference_sequence = ''.join(line.strip() for line in ref_lines if not line.startswith(">"))

if not reference_sequence:
    print("Warning: Reference sequence is empty.")
else:
    print(f"Reference sequence loaded successfully with length: {len(reference_sequence)}")

# Loop over VCF files like variants_r0.vcf, variants_r1.vcf, ...
for file_name in os.listdir(vcf_folder):
    if file_name.startswith("variants_") and file_name.endswith(".vcf"):
        vcf_path = os.path.join(vcf_folder, file_name)

        # Extract suffix (e.g., r0)
        suffix = file_name.replace("variants_", "").replace(".vcf", "")

        # Output FASTA file path
        output_fasta_path = os.path.join(output_folder, f"homo_ref_{suffix}.fasta")

        # Read sample names from the #CHROM line
        with open(vcf_path, "r") as vcf_file:
            for line in vcf_file:
                if line.startswith("#CHROM"):
                    columns = line.strip().split('\t')
                    sample_names = columns[9:]  # Columns after FORMAT
                    break
            else:
                print(f"Warning: No #CHROM line found in {vcf_path}")
                continue

        # Write to FASTA
        with open(output_fasta_path, "w") as output_fasta:
            for sample_name in sample_names:
                output_fasta.write(f">{sample_name}\n")
                output_fasta.write(f"{reference_sequence}\n")

        print(f"FASTA written: {output_fasta_path}")

print("✅ All FASTA files generated.")
