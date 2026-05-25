from Bio import SeqIO

# Mapping nucleotides to indices based on the given values
NUCLEOTIDE_TO_CONVERTED = {
    'A': '0', 'a': '0',
    'C': '5', 'c': '5',
    'G': 'a', 'g': 'a',
    'T': 'f', 't': 'f'
}

def convert_to_homozygous(sequence):
    """Convert the given sequence to the new format (A/a -> 0, C/c -> 5, G/g -> a, T/t -> f)."""
    homozygous_sequence = []
    for allele in sequence:
        if allele in NUCLEOTIDE_TO_CONVERTED:
            converted_nucleotide = NUCLEOTIDE_TO_CONVERTED[allele]  # Convert according to the new mapping
            homozygous_sequence.append(converted_nucleotide)
        else:
            raise ValueError(f"Unexpected allele: {allele}")
    return ''.join(homozygous_sequence)


def get_phased_genotype_index(parent1_index, parent2_index):
    """Convert phased genotype indices to a homozygous state using the formula."""
    if 0 <= parent1_index < 4 and 0 <= parent2_index < 4:
        return 4 * parent1_index + parent2_index
    else:
        raise RuntimeError("The parents should be canonical states.")

def process_fasta(gt16_fasta_path, homozygous_fasta_path):
    """Convert sequences in the FASTA file to the new homozygous format."""
    with open(gt16_fasta_path, "r") as gt16_file, open(homozygous_fasta_path, "w") as homozygous_file:
        for record in SeqIO.parse(gt16_file, "fasta"):
            # Convert the sequence to the new format
            original_sequence = str(record.seq)
            converted_sequence = convert_to_homozygous(original_sequence)

            # Write to the new FASTA file
            homozygous_file.write(f">{record.id}\n")
            homozygous_file.write(f"{converted_sequence}\n")

# File paths
ref_fasta = "/home/yil218/data/10K.fna"
homozygous_fasta = "/home/yil218/data/10K_homo.fna"

# Process the FASTA file and convert to the new format
process_fasta(ref_fasta, homozygous_fasta)

print(f"Converted GT16 FASTA to the new format: {homozygous_fasta}")
