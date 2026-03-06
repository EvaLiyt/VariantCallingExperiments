import os
import csv
import re
import statistics

# ---- Config ----
GENOTYPES = ["homo", "hetero"]
STATS_FOLDER = "./stats"
OUTPUT_FOLDER = "../summary"
OUTPUT_FILE = os.path.join(OUTPUT_FOLDER, "summary_method2.csv")


def parse_filename(filename):
    """Extract r, minD, and minVar from filename"""
    pattern = r"stats_r(\d+)_minD(\d+)_minVar(\d+)\.csv"
    match = re.match(pattern, filename)
    if match:
        return int(match.group(1)), int(match.group(2)), int(match.group(3))
    return None, None, None


def calculate_means(values):
    """Mean with empty list handling"""
    return statistics.mean(values) if values else 0.0


def read_overall_metrics(filepath):
    """Read overall metrics (Mean_TPR, Mean_FPR, ..., Mean_F1Score)"""
    metrics = {k: [] for k in [
        "TPR", "FPR", "FNR", "TNR",
        "TP", "FP", "TN", "FN",
        "Accuracy", "Precision", "F1Score"
    ]}

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
        return {f"Mean_{k}": 0.0 for k in metrics}

    return {f"Mean_{k}": calculate_means(v) for k, v in metrics.items()}


def read_all_columns(filepath):
    """
    Read *all* numeric columns from the stats file (e.g., 0_TP, 0_Precision, f_F1, etc.)
    and compute their means automatically.
    """
    means = {}
    with open(filepath, "r") as csvfile:
        reader = csv.DictReader(csvfile)
        headers = reader.fieldnames
        numeric_data = {h: [] for h in headers if h != "Cell"}

        for row in reader:
            for h in numeric_data:
                val = row[h]
                try:
                    numeric_data[h].append(float(val))
                except (ValueError, TypeError):
                    continue

    for h, vals in numeric_data.items():
        means[h] = calculate_means(vals) if vals else 0.0

    return means


def main():
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    if not os.path.exists(STATS_FOLDER):
        print(f"Stats folder {STATS_FOLDER} does not exist!")
        return

    summary_rows = []
    all_headers = set()

    for filename in os.listdir(STATS_FOLDER):
        if not filename.startswith("stats_r") or not filename.endswith(".csv"):
            continue

        r, minD, minVar = parse_filename(filename)
        if r is None:
            print(f"Could not parse filename: {filename}")
            continue

        filepath = os.path.join(STATS_FOLDER, filename)
        print(f"Processing {filename}...")

        # Compute overall mean metrics
        overall_means = read_overall_metrics(filepath)

        # Compute per-column means for everything in the file
        all_means = read_all_columns(filepath)

        # Merge both
        row = {
            "r_value": r,
            "minD_value": minD,
            "minVar_value": minVar,
        }
        row.update(overall_means)
        row.update(all_means)

        summary_rows.append(row)
        all_headers.update(row.keys())

    # --- Write summary CSV ---
    if not summary_rows:
        print("No valid stats files found.")
        return

    # Define column order
    fieldnames = ["r_value", "minD_value", "minVar_value"]
    mean_cols = [f"Mean_{k}" for k in [
        "TPR", "FPR", "FNR", "TNR", "TP", "FP", "TN", "FN", "Accuracy", "Precision", "F1Score"
    ]]
    other_cols = sorted(set(all_headers) - set(fieldnames) - set(mean_cols))
    fieldnames += mean_cols + other_cols

    with open(OUTPUT_FILE, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in sorted(summary_rows, key=lambda x: x["r_value"]):
            writer.writerow(row)

    print(f"\n Summary written to {OUTPUT_FILE}")
    print(f"Processed {len(summary_rows)} stats files successfully.")


if __name__ == "__main__":
    main()
