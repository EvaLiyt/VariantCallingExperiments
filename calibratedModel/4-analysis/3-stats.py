#!/usr/bin/env python3
import os
import csv
import re
from collections import defaultdict
from Bio import SeqIO

# --- Configuration ---
input_folder = "./comparison"
output_folder = "./stats"
site_folder = "./sites"
ref = "../2-variantCalling/merge"
ref_length = 10000

# Genotype mapping (degenerate codes -> canonical letter)
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
    gt = str(genotype)
    return genotype_map.get(gt, gt)


def read_genotype_from_files(r_folder):
    """Read genotype triplets from files in comparison/r_folder.
    Expects tab-separated lines: site <tab> ref_gen <tab> true_gen <tab> final_gen (or at least 4 cols).
    Returns dict keyed by (filename, site) -> (ref_gen, true_gen, final_gen) where genotypes are already mapped.
    """
    genotypes = {}
    print(f"[read_genotype_from_files] Reading genotypes from {r_folder}...")

    subfolder_path = os.path.join(input_folder, r_folder)
    if not os.path.isdir(subfolder_path):
        return genotypes

    for filename in sorted(os.listdir(subfolder_path)):
        file_path = os.path.join(subfolder_path, filename)
        if not os.path.isfile(file_path):
            continue

        # skip hidden/system files
        if filename.startswith("."):
            continue

        try:
            with open(file_path, 'r') as fh:
                for line in fh:
                    if not line.strip():
                        continue
                    if line.startswith("site"):
                        continue
                    cols = line.rstrip("\n").split('\t')
                    if len(cols) >= 4:
                        site = cols[0].strip()
                        # store mapped/upper-case genotypes
                        ref_genotype = map_genotype(cols[1].strip())
                        true_genotype = map_genotype(cols[2].strip())
                        final_genotype = map_genotype(cols[3].strip())
                        genotypes[(filename, site)] = (ref_genotype, true_genotype, final_genotype)
        except Exception as e:
            print(f"[read_genotype_from_files] Error reading {file_path}: {e}")

    return genotypes

def getRef(r_folder):
    """Load reference sequence for a given r_folder (expects ref/root_ref_{r_folder}.fasta).
    Returns a Seq object or an empty string on error.
    """
    ref_path = os.path.join(ref,f"root_ref_{r_folder}.fasta")
    ref_sequence = ""
    try:
        with open(ref_path) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                # if multiple sequences exist, keep last (mirrors original)
                ref_sequence = record.seq
    except Exception as e:
        print(f"[getRef] Error reading reference for {r_folder} at {ref_path}: {e}")
    return ref_sequence


def read_site_file(site_file_path):
    """Parse a site file with sections 'positive_sites' and 'negative_sites' followed by site lines.
    Returns two sets of site strings.
    """
    positive_sites = set()
    negative_sites = set()
    current_section = None

    try:
        with open(site_file_path, 'r') as site_file:
            for line in site_file:
                line = line.strip()
                if not line:
                    continue
                if line == "positive_sites":
                    current_section = "positive"
                elif line == "negative_sites":
                    current_section = "negative"
                elif current_section == "positive":
                    positive_sites.add(line)
                elif current_section == "negative":
                    negative_sites.add(line)
    except Exception as e:
        print(f"[read_site_file] Error reading site file {site_file_path}: {e}")

    return positive_sites, negative_sites

def compare_genotypes_with_sites(r_folder, genotypes):
    """For a given r_folder and genotypes dict, compute per-cell overall and per-genotype confusion stats."""
    print(f"[compare_genotypes_with_sites] Comparing genotypes with sites for {r_folder}...")
    results = {}

    subfolder_path = os.path.join(site_folder, r_folder)
    if not os.path.isdir(subfolder_path):
        return results

    # organize genotype data by cell_name -> {site: (ref, true, final)}
    genotype_by_cell = defaultdict(dict)
    for (filename, site), (ref_gen, true_gen, final_gen) in genotypes.items():
        cell_name = filename.split('_r')[0]
        genotype_by_cell[cell_name][site] = (ref_gen, true_gen, final_gen)

    # all sites as strings
    all_sites = {str(i) for i in range(1, ref_length + 1)}
    ref_sequence = getRef(r_folder)
	
    for filename in sorted(os.listdir(subfolder_path)):
        if filename.startswith("."):
            continue
        site_file_path = os.path.join(subfolder_path, filename)
        if not os.path.isfile(site_file_path):
            continue

        # parse positive/negative site lists
        positive_sites, negative_sites = read_site_file(site_file_path)
        cell_name = filename.split('_r')[0]

        stats = {'true_positive': 0, 'false_positive': 0, 'true_negative': 0, 'false_negative': 0}
        
        genotype_stats = {g: {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0} for g in {"homo", "hetero"}}
        cell_genotypes = genotype_by_cell.get(cell_name, {})
        homo = {"a", "0", "f", "5"}
        def is_homo(genotype):
            return genotype in homo
        
        for site in all_sites:
            if site in cell_genotypes:
                ref_genotype, true_genotype, final_genotype = cell_genotypes[site]
                # map the genotypes to keys we store counts under (avoid KeyError)
                true_key = map_genotype(true_genotype)
                final_key = map_genotype(final_genotype)
                is_match = true_key == final_key
                is_var = (ref_genotype is not None and true_key != ref_genotype)

                if site in positive_sites:
                    if is_var and is_match:
                        stats['true_positive'] += 1
                        if is_homo(true_key):
                            genotype_stats["homo"]['TP'] += 1
                        else:
                            genotype_stats["hetero"]['TP'] += 1
                    else:
                        stats['false_positive'] += 1
                        if is_homo(true_key):
                            genotype_stats["homo"]['FP'] += 1
                        else:
                            genotype_stats["hetero"]['FP'] += 1
                elif site in negative_sites:
                    if is_var:
                        stats['false_negative'] += 1
                        if is_homo(true_key):
                            genotype_stats["homo"]['FN'] += 1
                        else:
                            genotype_stats["hetero"]['FN'] += 1
                    else:
                        stats['true_negative'] += 1
                        if is_homo(true_key):
                            genotype_stats["homo"]['TN'] += 1
                        else:
                            genotype_stats["hetero"]['TN'] += 1
                else:
                    # unknown site -> treat as negative by default
                    if is_var:
                        stats['false_negative'] += 1
                        if is_homo(true_key):
                            genotype_stats["homo"]['FN'] += 1
                        else:
                            genotype_stats["hetero"]['FN'] += 1
                    else:
                        stats['true_negative'] += 1
                        if is_homo(true_key):
                            genotype_stats["homo"]['TN'] += 1
                        else:
                            genotype_stats["hetero"]['TN'] += 1
            else:
                stats['true_negative'] += 1
                
                idx = int(site) - 1
                genotype = str(ref_sequence[idx])
                mapped = map_genotype(genotype)
                if is_homo(mapped):
                    genotype_stats["homo"]['TN'] += 1
                else:
                    genotype_stats["hetero"]['TN'] += 1
                
        results[cell_name] = {'overall': stats, 'per_genotype': genotype_stats}
        
        # # check all genotypes measurements total number
        # totals = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}
        
        # # Sum across genotypes
        # for g, g_stats in genotype_stats.items():
        #     for key in totals:
        #         totals[key] += g_stats.get(key, 0)
        
        # print("\n--- Totals across all genotypes ---")
        # for k, v in totals.items():
        #     print(f"{k}: {v}")
    
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
    """Compute TPR, FPR, FNR, TNR, Precision, F1 per genotype"""
    results = {}
    for g, s in genotype_stats.items():
        tp, fp, tn, fn = s['TP'], s['FP'], s['TN'], s['FN']
        tpr = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        fpr = fp / (fp + tn) if (fp + tn) > 0 else 0.0
        fnr = fn / (tp + fn) if (tp + fn) > 0 else 0.0
        tnr = tn / (tn + fp) if (tn + fp) > 0 else 0.0
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        f1 = 2 * precision * tpr / (precision + tpr) if (precision + tpr) > 0 else 0.0
        results[g] = {'TP': tp, 'FP': fp, 'TN': tn, 'FN': fn,
                      'TPR': tpr, 'FPR': fpr, 'FNR': fnr, 'TNR': tnr,
                      'Precision': precision, 'F1': f1}
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
        print(f"[process_single_r] Error processing {r_folder}: {e}")
        return None

    if minD_value is None or minVar_value is None:
        print(f"[process_single_r] Could not extract minD/minVar for {r_folder}")
        return None

    try:
        genotypes = read_genotype_from_files(r_folder)
        if not genotypes:
            print(f"[process_single_r] No genotypes found for {r_folder}")
            return None

        cell_stats = compare_genotypes_with_sites(r_folder, genotypes)
        if not cell_stats:
            print(f"[process_single_r] No cell stats generated for {r_folder}")
            return None

        os.makedirs(output_folder, exist_ok=True)
        output_file = os.path.join(output_folder, f"stats_{r_folder}_minD{minD_value}_minVar{minVar_value}.csv")
        with open(output_file, 'w', newline='') as csvfile:
            base_fields = ['Cell', 'TP', 'FP', 'TN', 'FN', 'TPR', 'FPR', 'FNR', 'TNR', 'Accuracy', 'Precision', 'F1Score']

            genotype_fields = []
            for g in ["homo", "hetero"]:
                genotype_fields += [f"{g}_TP", f"{g}_FP", f"{g}_TN", f"{g}_FN",
                                    f"{g}_TPR", f"{g}_FPR", f"{g}_TNR", f"{g}_FNR",
                                    f"{g}_Precision", f"{g}_F1"]
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

        print(f"[process_single_r] Wrote {output_file}")
        return output_file

    except Exception as e:
        print(f"[process_single_r] Error processing r{r}: {e}")
        return None

def main():
    os.makedirs(output_folder, exist_ok=True)
    # clean previous outputs
    for filename in os.listdir(output_folder):
        file_path = os.path.join(output_folder, filename)
        if os.path.isfile(file_path):
            try:
                os.remove(file_path)
            except Exception:
                pass

    # detect available r folders under site_folder
    available_r_values = []
    for r in range(0, 100): 
        r_folder = f"r{r}"
        site_path = os.path.join(site_folder, r_folder)
        if os.path.exists(site_path):
            available_r_values.append(r)

    print(f"[main] Found {len(available_r_values)} r values to process")

    successful_files = []
    for r in available_r_values:
        result = process_single_r(r)
        if result:
            successful_files.append(result)

    print(f"[main] Successfully generated {len(successful_files)} stats files.")

if __name__ == "__main__":
    main()
