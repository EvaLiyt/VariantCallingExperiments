import os
import shutil
import re

# ======================
# PATH CONFIGURATION
# ======================
input_folder = "./files"
output_folder = "./comparison"
site_folder = "./sites"
params = "../../1-simulation/pileup"
vcf_folder = "../../2-variantCalling/filterMethods/method2/merge"

print("Current working directory:", os.getcwd())

# ======================
# PREPARE OUTPUT FOLDERS
# ======================
for folder in [output_folder, site_folder]:
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

# ======================
# GENOTYPE MAP
# ======================
genotype_map = {
    '1': 'M', '4': 'M',
    '2': 'R', '8': 'R',
    '3': 'W', 'c': 'W',
    '6': 'S', '9': 'S',
    '7': 'Y', 'd': 'Y',
    'b': 'K', 'e': 'K'
}

def map_genotype(genotype):
    if genotype is None:
        return genotype
    return genotype_map.get(str(genotype), genotype)

# ======================
# UTILS
# ======================
def clean_cell_name(sample_name):
    """Remove leading rX_ from sample name."""
    return re.sub(r"^r\d+_", "", sample_name)

def read_sequences_from_file(file_path):
    with open(file_path, 'r') as f:
        return [line.strip() for line in f.readlines()[1:]]

def compare_sequences(seq1, seq2, seq3):
    differing_sites = []
    for i in range(min(len(seq1), len(seq2), len(seq3))):
        site2 = map_genotype(seq2[i])
        site3 = map_genotype(seq3[i])
        if seq1[i] != seq2[i] or site2 != site3:
            differing_sites.append((i + 1, seq1[i], seq2[i], seq3[i]))
    return differing_sites

def parse_params(file_name):
    base = os.path.basename(file_name)
    if not base.startswith("parameters_"):
        return None, None, None

    r = base.replace("parameters_", "").replace(".txt", "")
    minD, minVar = None, None

    try:
        with open(file_name) as f:
            for line in f:
                if line.startswith("min_depth:"):
                    minD = line.split(":")[1].strip()
                elif line.startswith("min_var:"):
                    minVar = line.split(":")[1].strip()
    except FileNotFoundError:
        print(f"Missing parameter file: {file_name}")

    return r, minD, minVar

# ======================
# MAIN PER-r PROCESSING
# ======================
def process_r_folder(r_folder_name):

    sub_input_folder = os.path.join(input_folder, r_folder_name)
    sub_output_folder = os.path.join(output_folder, r_folder_name)
    sub_site_folder = os.path.join(site_folder, r_folder_name)

    os.makedirs(sub_output_folder, exist_ok=True)
    os.makedirs(sub_site_folder, exist_ok=True)

    params_file = os.path.join(params, r_folder_name, f"parameters_{r_folder_name}.txt")
    r, minD, minVar = parse_params(params_file)

    # ------------------
    # FASTA COMPARISON
    # ------------------
    if os.path.exists(sub_input_folder):
        for fasta_file in os.listdir(sub_input_folder):
            fasta_path = os.path.join(sub_input_folder, fasta_file)
            seqs = read_sequences_from_file(fasta_path)

            if len(seqs) != 3:
                continue

            diffs = compare_sequences(seqs[0], seqs[1], seqs[2])
            cell = fasta_file.split(".")[0]

            out_file = os.path.join(sub_output_folder, f"{cell}.tsv")
            with open(out_file, "w") as f:
                f.write("site\tref\ttrue\tfinal\n")
                for s in diffs:
                    f.write(f"{s[0]}\t{s[1]}\t{s[2]}\t{s[3]}\n")

    # ------------------
    # MERGED VCF (split by cell)
    # ------------------
    merged_vcf = os.path.join(vcf_folder, f"merged_{r_folder_name}.vcf")
    if not os.path.exists(merged_vcf):
        print(f"VCF not found: {merged_vcf}")
        return

    pos_sites = {}
    neg_sites = {}
    samples = []

    with open(merged_vcf) as vcf:
        gt_index = None

        for line in vcf:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                header = line.strip().split("\t")
                samples = header[9:]
                for s in samples:
                    pos_sites[s] = []
                    neg_sites[s] = []
                continue

            cols = line.strip().split("\t")
            fmt = cols[8].split(":")
            if "GT" not in fmt:
                continue
            gt_index = fmt.index("GT")
            pos = cols[1]

            for i, sample in enumerate(samples):
                fields = cols[9 + i].split(":")
                if gt_index >= len(fields):
                    continue
                gt = fields[gt_index]

                if gt == "./.":
                    continue
                elif gt in ("0/0", "0|0"):
                    neg_sites[sample].append(pos)
                else:
                    pos_sites[sample].append(pos)

    # ------------------
    # WRITE PER-CELL FILES
    # ------------------
    for sample in samples:
        cell = clean_cell_name(sample)
        cell = cell.replace("/", "_").replace(" ", "_")

        name = f"{cell}_r{r}_minD{minD}_minVar{minVar}.tsv"
        out_path = os.path.join(sub_site_folder, name)

        with open(out_path, "w") as f:
            f.write("positive_sites\n")
            f.write("\n".join(pos_sites[sample]) + "\n")
            f.write("negative_sites\n")
            f.write("\n".join(neg_sites[sample]) + "\n")

        print(f"Written: {out_path}")

# ======================
# ENTRY POINT
# ======================
def main():
    for i in range(100):
        process_r_folder(f"r{i}")

if __name__ == "__main__":
    main()
