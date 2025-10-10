import os
import csv
import re
import statistics

# ---- Config ----
GENOTYPES = ["0", "5", "K", "M", "R", "S", "W", "Y", "a", "f"]
STATS_FOLDER = "./stats"
OUTPUT_FOLDER = "./summary"
OUTPUT_FILE = os.path.join(OUTPUT_FOLDER, "summary.csv")


def parse_filename(filename):
    """Extract r, minD, and minVar from filename"""
    pattern = r"stats_r(\d+)_minD(\d+)_minVar(\d+)\.csv"
    match = re.match(pattern, filename)
    if match:
        r = int(match.group(1))
        minD = int(match.group(2))
        minVar = int(match.group(3))
        return r, minD, minVar
    else:
        return None, None, None


def read_genotype_rates(filepath):
    """
    Read per-genotype rate columns (*_TPR, *_FPR, *_TNR, *_FNR)
    and compute average across all cells.
    """
    genotype_data = {g: {"TPR": [], "FPR": [], "TNR": [], "FNR": []} for g in GENOTYPES}

    try:
        with open(filepath, "r") as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                for g in GENOTYPES:
                    for metric in ["TPR", "FPR", "TNR", "FNR"]:
                        col = f"{g}_{metric}"
                        if col in row and row[col] not in [None, ""]:
                            try:
                                genotype_data[g][metric].append(float(row[col]))
                            except ValueError:
                                continue
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return None

    # Compute means
    genotype_means = {}
    for g in GENOTYPES:
        for metric in ["TPR", "FPR", "TNR", "FNR"]:
            values = genotype_data[g][metric]
            mean_val = statistics.mean(values) if values else 0.0
            genotype_means[f"{g}_{metric}"] = mean_val

    return genotype_means


def calculate_means(values):
    """Mean with empty list handling"""
    if not values:
        return 0.0
    return statistics.mean(values)


def read_overall_metrics(filepath):
    """Read overall TPR/FPR/FNR/TNR/Accuracy/Precision/F1Score"""
    metrics = {k: [] for k in ["TPR", "FPR", "FNR", "TNR", "Accuracy", "Precision", "F1Score"]}

    try:
        with open(filepath, "r") as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                for key in metrics.keys():
                    try:
                        metrics[key].append(float(row[key]))
                    except (ValueError, KeyError):
                        continue
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return {k: [] for k in metrics}

    mean_metrics = {f"Mean_{k}": calculate_means(v) for k, v in metrics.items()}
    return mean_metrics


def main():
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    if not os.path.exists(STATS_FOLDER):
        print(f"Stats folder {STATS_FOLDER} does not exist!")
        return

    summary_rows = []

    for filename in os.listdir(STATS_FOLDER):
        if not filename.startswith("stats_r") or not filename.endswith(".csv"):
            continue

        r, minD, minVar = parse_filename(filename)
        if r is None:
            print(f"Could not parse filename: {filename}")
            continue

        filepath = os.path.join(STATS_FOLDER, filename)
        print(f"Processing {filename}...")

        # Read overall metrics
        overall_means = read_overall_metrics(filepath)

        # Read genotype-specific averages
        genotype_means = read_genotype_rates(filepath)
        if not genotype_means:
            print(f"No genotype metrics found in {filename}")
            continue

        row = {
            "r_value": r,
            "minD_value": minD,
            "minVar_value": minVar,
        }
        row.update(overall_means)
        row.update(genotype_means)

        summary_rows.append(row)

    # Sort results by r
    summary_rows.sort(key=lambda x: x["r_value"])

    # Prepare header
    if summary_rows:
        fieldnames = ["r_value", "minD_value", "minVar_value"]
        fieldnames += [f"Mean_{k}" for k in ["TPR", "FPR", "FNR", "TNR", "Accuracy", "Precision", "F1Score"]]
        for g in GENOTYPES:
            fieldnames += [f"{g}_{m}" for m in ["TPR", "FPR", "TNR", "FNR"]]

        # Write output CSV
        with open(OUTPUT_FILE, "w", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for row in summary_rows:
                writer.writerow(row)

        print(f"\n✅ Summary written to {OUTPUT_FILE}")
        print(f"Processed {len(summary_rows)} stats files successfully.")
    else:
        print("No valid stats files found.")


if __name__ == "__main__":
    main()
