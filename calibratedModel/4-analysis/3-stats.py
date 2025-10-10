import os
import shutil
import csv
import re
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import multiprocessing

# Define folder paths
input_folder = "./comparison"
output_folder = "./stats"
site_folder = "./sites"
ref_length = 10000

# Genotype mapping
genotype_map = {
    '1': 'M', '4': 'M',
    '2': 'R', '8': 'R',
    '3': 'W', 'c': 'W',
    '6': 'S', '9': 'S',
    '7': 'Y', 'd': 'Y',
    'b': 'K', 'e': 'K'
}

# Include unmapped (non-degenerate) genotypes as themselves
all_genotypes = sorted(set(genotype_map.values()) | {'0', '5', 'a', 'f'})


def map_genotype(genotype):
    return genotype_map.get(genotype, genotype)


def read_genotype_from_files(r_folder):
    """Optimized genotype reading with batch processing"""
    genotypes = {}
    print(f"Reading genotypes from {r_folder}...")

    subfolder_path = os.path.join(input_folder, r_folder)
    if not os.path.isdir(subfolder_path):
        return genotypes

    for filename in os.listdir(subfolder_path):
        file_path = os.path.join(subfolder_path, filename)
        if not os.path.isfile(file_path):
            continue

        try:
            with open(file_path, 'r') as file:
                lines = file.readlines()
                for line in lines:
                    if line.startswith("site"):
                        continue

                    cols = line.strip().split('\t')
                    if len(cols) >= 4:
                        site = cols[0]
                        ref_genotype = map_genotype(cols[1])
                        true_genotype = map_genotype(cols[2])
                        final_genotype = map_genotype(cols[3])
                        genotypes[(filename, site)] = (ref_genotype, true_genotype, final_genotype)
        except Exception as e:
            print(f"Error reading {file_path}: {e}")

    return genotypes


def read_site_file(site_file_path):
    positive_sites = set()
    negative_sites = set()
    current_section = None

    try:
        with open(site_file_path, 'r') as site_file:
            for line in site_file:
                line = line.strip()
                if line == "positive_sites":
                    current_section = "positive"
                elif line == "negative_sites":
                    current_section = "negative"
                elif current_section == "positive" and line:
                    positive_sites.add(line)
                elif current_section == "negative" and line:
                    negative_sites.add(line)
    except Exception as e:
        print(f"Error reading site file {site_file_path}: {e}")

    return positive_sites, negative_sites


def compare_genotypes_with_sites(r_folder, genotypes):
    """Optimized genotype comparison with per-genotype tracking"""
    print(f"Comparing genotypes with sites for {r_folder}...")
    results = {}

    subfolder_path = os.path.join(site_folder, r_folder)
    if not os.path.isdir(subfolder_path):
        return results

    genotype_by_cell = defaultdict(dict)
    for (filename, site), (ref_gen, true_gen, final_gen) in genotypes.items():
        cell_name = filename.split('_r')[0]
        genotype_by_cell[cell_name][site] = (ref_gen, true_gen, final_gen)

    all_sites = {str(i) for i in range(1, ref_length + 1)}

    for filename in os.listdir(subfolder_path):
        site_file_path = os.path.join(subfolder_path, filename)
        if not os.path.isfile(site_file_path):
            continue

        positive_sites, negative_sites = read_site_file(site_file_path)
        cell_name = filename.split('_r')[0]

        stats = {'true_positive': 0, 'false_positive': 0, 'true_negative': 0, 'false_negative': 0}
        # NEW: per-genotype confusion stats
        genotype_stats = {g: {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0} for g in all_genotypes}

        cell_genotypes = genotype_by_cell.get(cell_name, {})

        for site in all_sites:
            if site in cell_genotypes:
                ref_genotype, true_genotype, final_genotype = cell_genotypes[site]
                is_match = true_genotype == final_genotype
                is_var = ref_genotype != true_genotype

                # Overall metrics
                if site in positive_sites:
                    if is_var and is_match:
                        stats['true_positive'] += 1
                    else:
                        stats['false_positive'] += 1
                elif site in negative_sites:
                    if is_var:
                        stats['false_negative'] += 1
                    else:
                        stats['true_negative'] += 1
                else:
                    if is_var:
                        stats['false_negative'] += 1
                    else:
                        stats['true_negative'] += 1

                # --- Per-genotype metrics ---
                if true_genotype in all_genotypes:
                    if final_genotype == true_genotype:
                        genotype_stats[true_genotype]['TP'] += 1
                    else:
                        genotype_stats[true_genotype]['FN'] += 1

                if final_genotype in all_genotypes and final_genotype != true_genotype:
                    genotype_stats[final_genotype]['FP'] += 1

                for g in all_genotypes:
                    if g not in (true_genotype, final_genotype):
                        genotype_stats[g]['TN'] += 1
            else:
                # No genotype data
                stats['true_negative'] += 1
                for g in all_genotypes:
                    genotype_stats[g]['TN'] += 1

        results[cell_name] = {'overall': stats, 'per_genotype': genotype_stats}

    return results


def calc_metrics(stats):
    tp, fp, tn, fn = stats['true_positive'], stats['false_positive'], stats['true_negative'], stats['false_negative']
    tpr = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    fpr = fp / (fp + tn) if (fp + tn) > 0 else 0.0
    fnr = fn / (tp + fn) if (tp + fn) > 0 else 0.0
    tnr = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    accuracy = (tp + tn) / (tp + tn + fp + fn) if (tp + tn + fp + fn) > 0 else 0.0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    f1 = 2 * precision * tpr / (tpr + precision) if (tpr + precision) > 0 else 0.0
    return tpr, fpr, fnr, tnr, accuracy, precision, f1


def calc_per_genotype_metrics(genotype_stats):
    """Compute TPR, Precision, F1 per genotype"""
    results = {}
    for g, s in genotype_stats.items():
        tp, fp, tn, fn = s['TP'], s['FP'], s['TN'], s['FN']
        tpr = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        fpr = fp / (fp + tn) if (fp + tn) > 0 else 0.0
        fnr = fn / (tp + fn) if (tp + fn) > 0 else 0.0
        tnr = tn / (tn + fp) if (tn + fp) > 0 else 0.0
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        f1 = 2 * precision * tpr / (precision + tpr) if (precision + tpr) > 0 else 0.0
        results[g] = {'TP': tp, 'FP': fp, 'TN': tn, 'FN': fn, 'TPR': tpr, 'FPR':fpr, 'FNR':fnr, 'TNR':tnr,'Precision': precision, 'F1': f1}
    return results


def sorted_cell_names(cell_names):
    number_pattern = re.compile(r'\d+')

    def extract_number(name):
        match = number_pattern.search(name)
        return int(match.group()) if match else 0

    clade = [name for name in cell_names if "clade" in name]
    normal = [name for name in cell_names if "clade" not in name]

    return sorted(normal, key=extract_number) + sorted(clade, key=extract_number)


def process_single_r(r):
    r_folder = f"r{r}"
    site_path = os.path.join(site_folder, r_folder)
    if not os.path.exists(site_path):
        return None

    minD_value = minVar_value = None
    try:
        for filename in os.listdir(site_path):
            match = re.search(r"minD(\d+)_minVar(\d+)", filename)
            if match:
                minD_value = int(match.group(1))
                minVar_value = int(match.group(2))
                break
    except Exception as e:
        print(f"Error processing {r_folder}: {e}")
        return None

    if minD_value is None or minVar_value is None:
        print(f"Could not extract minD/minVar for r{r}")
        return None

    try:
        genotypes = read_genotype_from_files(r_folder)
        if not genotypes:
            print(f"No genotypes found for {r_folder}")
            return None

        cell_stats = compare_genotypes_with_sites(r_folder, genotypes)
        if not cell_stats:
            print(f"No cell stats generated for {r_folder}")
            return None

        output_file = os.path.join(output_folder, f"stats_{r_folder}_minD{minD_value}_minVar{minVar_value}.csv")
        with open(output_file, 'w', newline='') as csvfile:
            # --- Define all field names ---
            base_fields = ['Cell', 'TP', 'FP', 'TN', 'FN', 'TPR', 'FPR', 'FNR', 'TNR', 'Accuracy', 'Precision', 'F1Score']
            genotype_fields = []
            for g in all_genotypes:
                genotype_fields += [f"{g}_TP", f"{g}_FP", f"{g}_TN", f"{g}_FN", f"{g}_TPR", f"{g}_FPR",f"{g}_TNR",f"{g}_FNR", f"{g}_Precision", f"{g}_F1"]
            fieldnames = base_fields + genotype_fields

            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()

            for cell in sorted_cell_names(cell_stats.keys()):
                stats = cell_stats[cell]['overall']
                tpr, fpr, fnr, tnr, accuracy, precision, f1 = calc_metrics(stats)
                row = {
                    'Cell': cell,
                    'TP': stats['true_positive'],
                    'FP': stats['false_positive'],
                    'TN': stats['true_negative'],
                    'FN': stats['false_negative'],
                    'TPR': tpr,
                    'FPR': fpr,
                    'FNR': fnr,
                    'TNR': tnr,
                    'Accuracy': accuracy,
                    'Precision': precision,
                    'F1Score': f1
                }

                per_geno = calc_per_genotype_metrics(cell_stats[cell]['per_genotype'])
                for g, m in per_geno.items():
                    row[f"{g}_TP"] = m['TP']
                    row[f"{g}_FP"] = m['FP']
                    row[f"{g}_TN"] = m['TN']
                    row[f"{g}_FN"] = m['FN']
                    row[f"{g}_TPR"] = m['TPR']
                    row[f"{g}_FPR"] = m['FPR']
                    row[f"{g}_TNR"] = m['TNR']
                    row[f"{g}_FNR"] = m['FNR']
                    row[f"{g}_Precision"] = m['Precision']
                    row[f"{g}_F1"] = m['F1']

                writer.writerow(row)

        print(f"Wrote {output_file}")
        return output_file

    except Exception as e:
        print(f"Error processing r{r}: {e}")
        return None


def main():
    os.makedirs(output_folder, exist_ok=True)
    for filename in os.listdir(output_folder):
        file_path = os.path.join(output_folder, filename)
        if os.path.isfile(file_path):
            os.remove(file_path)

    available_r_values = []
    for r in range(0, 100):
        r_folder = f"r{r}"
        site_path = os.path.join(site_folder, r_folder)
        if os.path.exists(site_path):
            available_r_values.append(r)

    print(f"Found {len(available_r_values)} r values to process")

    successful_files = []
    for r in available_r_values:
        result = process_single_r(r)
        if result:
            successful_files.append(result)

    print(f"Successfully generated {len(successful_files)} stats files.")


if __name__ == "__main__":
    main()
