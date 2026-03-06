#!/usr/bin/env python3
import os
import re
import csv
from collections import defaultdict

vcf_directory = "../../2-variantCalling/filterMethods/method3/merge"
param_directory = "../../1-simulation/pileup"
output_directory = "../summary"
output_csv = os.path.join(output_directory, "var_counts_method3.csv")

def read_parameters(r_value):
    path = os.path.join(param_directory, f"r{r_value}", f"parameters_r{r_value}.txt")
    if not os.path.exists(path):
        return None, None
    with open(path) as f:
        lines = [l.strip() for l in f if l.strip()]
    return (
        lines[0].split(":")[1].strip(),
        lines[1].split(":")[1].strip(),
    )

def classify_sample(sample_name):
    # strip rX_
    name = re.sub(r"^r\d+_", "", sample_name)
    return "tumour" if name.startswith("clade_") else "normal"

def summarize_merged_vcf(filepath):
    samples = None
    counts = None

    with open(filepath) as f:
        for line in f:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                header = line.rstrip().split("\t")
                samples = header[9:]
                counts = {s: {"total": 0, "zero_zero": 0, "non_zero": 0} for s in samples}
                continue

            if line.startswith("#"):
                continue

            cols = line.rstrip().split("\t")
            fmt = cols[8].split(":")
            if "GT" not in fmt:
                continue
            gt_idx = fmt.index("GT")

            for s, sf in zip(samples, cols[9:]):
                gt = sf.split(":")[gt_idx]

                if gt in ("./.", ".|."):
                    continue

                counts[s]["total"] += 1
                if gt in ("0/0", "0|0"):
                    counts[s]["zero_zero"] += 1
                else:
                    counts[s]["non_zero"] += 1

    return samples, counts

variant_stats = defaultdict(lambda: {
    "normal_total_variants": [],
    "normal_zero_zero": [],
    "normal_non_zero": [],
    "tumour_total_variants": [],
    "tumour_zero_zero": [],
    "tumour_non_zero": []
})

for fname in sorted(os.listdir(vcf_directory)):
    m = re.match(r"merged_r(\d+)\.vcf$", fname)
    if not m:
        continue

    r = m.group(1)
    minD, minVar = read_parameters(r)
    if minD is None:
        continue

    samples, per_sample = summarize_merged_vcf(
        os.path.join(vcf_directory, fname)
    )

    key = (r, minD, minVar)
    stats = variant_stats[key]

    for s in samples:
        cls = classify_sample(s)
        if cls == "tumour":
            stats["tumour_total_variants"].append(per_sample[s]["total"])
            stats["tumour_zero_zero"].append(per_sample[s]["zero_zero"])
            stats["tumour_non_zero"].append(per_sample[s]["non_zero"])
        else:
            stats["normal_total_variants"].append(per_sample[s]["total"])
            stats["normal_zero_zero"].append(per_sample[s]["zero_zero"])
            stats["normal_non_zero"].append(per_sample[s]["non_zero"])

    print(f"merged_r{r}.vcf: samples={len(samples)}")

os.makedirs(output_directory, exist_ok=True)

def avg(x):
    return round(sum(x) / len(x), 2) if x else 0.0

with open(output_csv, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=[
        "r_value", "minD_value", "minVar_value",
        "avg_total_variants_normal", "avg_zero_zero_normal", "avg_non_zero_normal",
        "avg_total_variants_tumour", "avg_zero_zero_tumour", "avg_non_zero_tumour"
    ])
    writer.writeheader()

    for (r, minD, minVar), s in sorted(variant_stats.items(), key=lambda x: int(x[0][0])):
        writer.writerow({
            "r_value": r,
            "minD_value": minD,
            "minVar_value": minVar,
            "avg_total_variants_normal": avg(s["normal_total_variants"]),
            "avg_zero_zero_normal": avg(s["normal_zero_zero"]),
            "avg_non_zero_normal": avg(s["normal_non_zero"]),
            "avg_total_variants_tumour": avg(s["tumour_total_variants"]),
            "avg_zero_zero_tumour": avg(s["tumour_zero_zero"]),
            "avg_non_zero_tumour": avg(s["tumour_non_zero"]),
        })

print("Done.")
