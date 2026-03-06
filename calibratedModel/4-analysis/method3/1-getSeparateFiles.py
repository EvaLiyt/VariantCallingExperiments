import os
import shutil
from Bio import SeqIO

output_folder = "./files"
# Check if the folder exists, if not create it
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
    print(f"Directory '{output_folder}' created.")
else:
    # If it exists, delete all files in it
    for filename in os.listdir(output_folder):
        file_path = os.path.join(output_folder, filename)
        try:
            if os.path.isfile(file_path):
                os.remove(file_path)
                print(f"Removed file: {file_path}")
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
                print(f"Removed directory: {file_path}")
        except Exception as e:
            print(f"Error removing {file_path}: {e}")

def read_nexus_sequences(nexus_file):
    """Reads sequences from a Nexus file."""
    with open(nexus_file, 'r') as f:
        sequences = {}
        in_matrix = False
        for line in f:
            line = line.strip()
            if line.lower() == "matrix":  # Start of the sequence block
                in_matrix = True
                continue
            elif line == ";":  # End of the sequence block
                in_matrix = False
            elif in_matrix:
                # Split line into taxon name and sequence
                parts = line.split()
                if len(parts) == 2:
                    taxon, seq = parts
                    sequences[taxon] = seq
        return sequences

def read_fasta_sequences(fasta_file):
    """Reads sequences from a FASTA file."""
    sequences = {}
    try:
        sequences = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
    return sequences

def main(nexus_file, fasta_file, reference_nexus, r_folder):
    # Create a subfolder for each 'r' inside the output folder
    output_r_folder = os.path.join(output_folder, r_folder)
    if not os.path.exists(output_r_folder):
        os.makedirs(output_r_folder)
        print(f"Directory '{output_r_folder}' created.")

    # Read sequences from the main Nexus and FASTA files
    nexus_sequences = read_nexus_sequences(nexus_file)
    fasta_sequences = read_fasta_sequences(fasta_file)

    # Read reference sequence from the reference Nexus file
    try:
        reference_sequences = read_nexus_sequences(reference_nexus)
        # You can change 'RefSeq' to the exact sequence name in your reference nexus file
        ref_key = list(reference_sequences.keys())[0]  # Get the first (or specify)
        reference_seq = reference_sequences[ref_key]
        #print(f"Using reference sequence '{ref_key}': {reference_seq}")
    except Exception as e:
        print(f"Error reading reference Nexus file: {e}")
        return

    # Generate output files
    for seq_name in nexus_sequences.keys():
        output_file_name = f"{seq_name}_{r_folder}.fasta"
        output_file_path = os.path.join(output_r_folder, output_file_name)
        print(f"Generating file: {output_file_path}")
        try:
            with open(output_file_path, 'w') as output_file:
                output_file.write(f">{seq_name}\n")
                output_file.write(f"{reference_seq}\n")
                output_file.write(f"{nexus_sequences[seq_name]}\n")
                output_file.write(f"{fasta_sequences.get(seq_name, 'N/A')}\n")
        except Exception as e:
            print(f"Error writing to file {output_file_path}: {e}")

if __name__ == "__main__":
    for i in range(100):  # Loop through r0 to r99
        r_folder = f"r{i}"
        nexus_file = f"../../1-simulation/nexus/sequenceGeneration_r{i}_A.nexus"
        fasta_file = f"../../2-variantCalling/filterMethods/method3/alignment/alignment_r{i}.fasta"
        reference_nexus = f"../../1-simulation/root/sequenceGeneration_r{i}_D.nexus"
        main(nexus_file, fasta_file, reference_nexus, r_folder)
