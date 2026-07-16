from Bio import SeqIO
import shutil
import os

input_folder = './merge'
output_folder = './alignment'

if os.path.exists(output_folder):
    print(f"Deleting existing folder: {output_folder}")
    shutil.rmtree(output_folder)
os.makedirs(output_folder)
print(f"Created the output folder: {output_folder}")

ref_files = [f for f in os.listdir(input_folder) if f.startswith("root_ref") and f.endswith(".fasta")]

for ref_fasta in ref_files:
    suffix = '_'.join(ref_fasta.split('_')[2:]).split('.')[0]
    ref_fasta_path = os.path.join(input_folder, ref_fasta)
    input_fasta_path = os.path.join(input_folder, f"raw_alignment_{suffix}.fasta")
    output_fasta_path = os.path.join(output_folder, f"alignment_{suffix}.fasta")

    if not os.path.exists(input_fasta_path):
        print(f"Warning: {input_fasta_path} not found.")
        continue

    # Load reference sequences into memory
    ref_dict = {rec.id: str(rec.seq) for rec in SeqIO.parse(ref_fasta_path, "fasta")}

    # Process alignment file
    modified_records = []
    for rec in SeqIO.parse(input_fasta_path, "fasta"):
        ref_seq = ref_dict.get(rec.id)
        if ref_seq:
            # Replace all ? in a list-based approach (fast)
            seq_list = list(str(rec.seq))
            for i, base in enumerate(seq_list):
                if base == '?':
                    if i < len(ref_seq):
                        seq_list[i] = ref_seq[i]
            new_seq = ''.join(seq_list)
        else:
            new_seq = str(rec.seq)  # No replacement if reference missing

        # Rename ID
        if 'clade' in rec.id:
            rec.id = '_'.join(rec.id.split('_')[:2])
        else:
            rec.id = rec.id.split('_')[0]
        rec.description = ""  # Clear description
        rec.seq = new_seq  # Store string directly

        modified_records.append(rec)

    # Write all at once
    with open(output_fasta_path, "w") as out:
        for rec in modified_records:
            out.write(f">{rec.id}\n{rec.seq}\n")

    print(f"Saved: {output_fasta_path}")

print("Done.")
