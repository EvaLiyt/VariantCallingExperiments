import os
import shutil
import csv
from Bio import SeqIO
import re

# Define folder paths
input_folder = "./comparison"
output_folder = "./stats"
site_folder = "./sites"
vcf_folder = "../2-variantCalling/filtered"
ref_folder = "../1-simulation/root"
ref = "/home/yil218/data/10K.fna"
ref_length = len(next(SeqIO.parse(ref, "fasta")).seq)

# Ensure output folder exists and clear it
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
else:
    for filename in os.listdir(output_folder):
        file_path = os.path.join(output_folder, filename)
        try:
            if os.path.isfile(file_path):
                os.remove(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print(f"Error removing {file_path}: {e}")

# Genotype mapping
genotype_map = {
    '1': 'M', '4': 'M',
    '2': 'R', '8': 'R',
    '3': 'W', 'c': 'W',
    '6': 'S', '9': 'S',
    '7': 'Y', 'd': 'Y',
    'b': 'K', 'e': 'K'
}


# Function to apply genotype mapping
def map_genotype(genotype):
    return genotype_map.get(genotype, genotype)


# Function to read genotypes from input files
def read_genotype_from_files():
    genotypes = {}
    print("Reading genotypes from files in the comparison folder...")
    for subfolder in os.listdir(input_folder):
        subfolder_path = os.path.join(input_folder, subfolder)
        if os.path.isdir(subfolder_path):
            print(f"Reading from subfolder: {subfolder}")
            for filename in os.listdir(subfolder_path):
                if filename.endswith("rx_pairwise.tsv"):
                    print(f"Skipping file: {filename}")
                    continue
                file_path = os.path.join(subfolder_path, filename)
                if os.path.isfile(file_path):
                    print(f"Processing file: {filename}")
                    with open(file_path, 'r') as file:
                        for line in file.readlines():
                            if line.startswith("site"):
                                continue
                            cols = line.strip().split('\t')
                            if len(cols) < 4:
                                print(f"Warning: Skipping malformed line in {filename}: {line.strip()}")
                                continue
                            site = cols[0]
                            true_genotype = map_genotype(cols[2])
                            final_genotype = map_genotype(cols[3])
                            genotypes[(filename, site)] = (true_genotype, final_genotype)
    return genotypes


def compare_genotypes_with_sites(genotypes):
    results = {}
    print("Comparing genotypes with sites...")

    # Cache ref lengths to avoid repeated loading
    ref_lengths = {}

    for subfolder in os.listdir(site_folder):
        subfolder_path = os.path.join(site_folder, subfolder)
        if os.path.isdir(subfolder_path):
            print(f"  Processing subfolder: {subfolder}")
            for filename in os.listdir(subfolder_path):
                site_file_path = os.path.join(subfolder_path, filename)
                if os.path.isfile(site_file_path):
                    print(f"Reading site file: {filename}")
                    positive_sites = set()
                    negative_sites = set()
                    current_section = None
                    with open(site_file_path, 'r') as site_file:
                        for line in site_file:
                            line = line.strip()
                            if line == "positive_sites":
                                current_section = "positive"
                                continue
                            elif line == "negative_sites":
                                current_section = "negative"
                                continue
                            if current_section == "positive":
                                positive_sites.add(line)
                            elif current_section == "negative":
                                negative_sites.add(line)

                    match = re.search(r"r(\d+)_minD(\d+)_minVar(\d+)", filename)
                    if match:
                        r_value = int(match.group(1))
                        minD_value = int(match.group(2))
                        minVar_value = int(match.group(3))
                    else:
                        print(f"Filename {filename} does not match expected pattern.")
                        continue

                    cell_name = filename.split('_r')[0]

                    cell_genotypes = {
                        site: (true_genotype, final_genotype)
                        for (file, site), (true_genotype, final_genotype) in genotypes.items()
                        if file.split("_r")[0] == cell_name
                    }

                    stats = {'true_positive': 0, 'false_positive': 0, 'true_negative': 0, 'false_negative': 0}
                    details = []

                    comparison_sites = set(cell_genotypes.keys())
                    all_sites = positive_sites | negative_sites
                    missing_in_comparison = all_sites - comparison_sites

                    for site in missing_in_comparison:
                        stats['true_negative'] += 1
                        details.append((site, "NA", "true_negative"))

                    for site, (true_genotype, final_genotype) in cell_genotypes.items():
                        mapped_true_genotype = map_genotype(true_genotype)
                        mapped_final_genotype = map_genotype(final_genotype)

                        if site in positive_sites:
                            if mapped_final_genotype == mapped_true_genotype:
                                stats['true_positive'] += 1
                            else:
                                stats['false_positive'] += 1
                                details.append((site, final_genotype, "false_positive"))
                        elif site in negative_sites:
                            if mapped_final_genotype == mapped_true_genotype:
                                stats['true_negative'] += 1
                            else:
                                stats['false_negative'] += 1
                                details.append((site, true_genotype, "false_negative"))
                        else:
                            stats['false_negative'] += 1

                    considered_num = stats['true_positive'] + stats['false_positive'] + stats['true_negative'] + stats['false_negative']
                    stats['true_negative'] = ref_length - considered_num

                    results[(r_value, minD_value, minVar_value, cell_name)] = (stats, details)

    return results


def calc_true_positive_rate(results):
    tpr_results = {}
    for key, (stats, _) in results.items():
        tp = stats['true_positive']
        fn = stats['false_negative']
        tpr = tp / (tp + fn) if (tp + fn) > 0 else None
        tpr_results[key] = tpr
    return tpr_results


def calc_false_positive_rate(results):
    fpr_results = {}
    for key, (stats, _) in results.items():
        fp = stats['false_positive']
        tn = stats['true_negative']
        fpr = fp / (fp + tn) if (fp + tn) > 0 else None
        fpr_results[key] = fpr
    return fpr_results


def calc_precision(results):
    precision_results = {}
    for key, (stats, _) in results.items():
        tp = stats['true_positive']
        fp = stats['false_positive']
        precision = tp / (tp + fp) if (tp + fp) > 0 else None
        precision_results[key] = precision
    return precision_results


def sorted_cell_names(cell_names):
    clade_names = sorted([name for name in cell_names if "clade" in name],
                         key=lambda x: int(re.search(r'\d+', x).group()))
    normal_names = sorted([name for name in cell_names if "clade" not in name],
                          key=lambda x: int(re.search(r'\d+', x).group()))
    return normal_names + clade_names


def main():
    print("Starting the analysis...")

    genotypes = read_genotype_from_files()
    print("Genotypes loaded, now comparing with sites...")

    results = compare_genotypes_with_sites(genotypes)
    print("Comparison complete. Results generated.")

    tpr_results = calc_true_positive_rate(results)
    fpr_results = calc_false_positive_rate(results)
    precision_results = calc_precision(results)

    summary_path = os.path.join(output_folder, "summary.csv")
    summary_exists = os.path.isfile(summary_path)

    with open(summary_path, 'a', newline='') as summary_file:
        summary_writer = csv.writer(summary_file)
        if not summary_exists:
            summary_writer.writerow(['r_value', 'minD_value', 'minVar_value', 'Mean_TPR', 'Mean_FPR', 'Mean_Precision'])

        for (r_value, minD_value, minVar_value) in sorted({(r, d, v) for r, d, v, _ in results.keys()}):
            filtered_results = [
                (tpr_results[k], fpr_results[k], precision_results[k])
                for k in results.keys() if k[0] == r_value and k[1] == minD_value and k[2] == minVar_value
            ]

            mean_tpr = sum(t[0] for t in filtered_results if t[0] is not None) / len(filtered_results)
            mean_fpr = sum(t[1] for t in filtered_results if t[1] is not None) / len(filtered_results)
            mean_precision = sum(t[2] for t in filtered_results if t[2] is not None) / len(filtered_results)

            summary_writer.writerow([r_value, minD_value, minVar_value, mean_tpr, mean_fpr, mean_precision])
            print(f"Summary for r={r_value}, minD={minD_value}, minVar={minVar_value} written.")

    print("Writing cell-specific files...")
    for (r_value, minD_value, minVar_value) in sorted({(r, d, v) for r, d, v, _ in results.keys()}):
        csv_filename = f"stats_r{r_value}_minD{minD_value}_minVar{minVar_value}.csv"
        csv_path = os.path.join(output_folder, csv_filename)

        with open(csv_path, 'w', newline='') as csvfile:
            fieldnames = ['Cell', 'TP', 'FP', 'TN', 'FN', 'TPR', 'FPR', 'Precision']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()

            cell_names = [cell_name for (r, d, v, cell_name) in results.keys()
                          if r == r_value and d == minD_value and v == minVar_value]
            sorted_names = sorted_cell_names(cell_names)

            for cell_name in sorted_names:
                (stats, _) = results[(r_value, minD_value, minVar_value, cell_name)]
                writer.writerow({
                    'Cell': cell_name,
                    'TP': stats['true_positive'],
                    'FP': stats['false_positive'],
                    'TN': stats['true_negative'],
                    'FN': stats['false_negative'],
                    'TPR': tpr_results[(r_value, minD_value, minVar_value, cell_name)],
                    'FPR': fpr_results[(r_value, minD_value, minVar_value, cell_name)],
                    'Precision': precision_results[(r_value, minD_value, minVar_value, cell_name)],
                })
        print(f"File {csv_filename} written.")

    print("All files have been successfully generated!")


if __name__ == "__main__":
    main()
