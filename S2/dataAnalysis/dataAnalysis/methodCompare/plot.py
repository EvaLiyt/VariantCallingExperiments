#!/usr/bin/env python3
"""
plot_results_enhanced.py
────────────────────────
Enhanced scatter plotting from per_replicate.csv with:
  ✓ Individual replicates as light dots (scattered around x-position)
  ✓ Darker, larger dots showing the mean at each D/V
  ✓ Error bars (±1 std) for each D/V group
  ✓ Clean, publication-ready styling
  ✓ All 4 metrics × 3 grouping strategies
"""

import os
import sys
import csv
from collections import defaultdict
import numpy as np

try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.patches import Rectangle
except ImportError:
    print("[error] matplotlib not installed. Run: pip install matplotlib --break-system-packages")
    sys.exit(1)

# ── Config ─────────────────────────────────────────────────────────────────────
CSV_PATH = "./results/per_replicate.csv"
OUTPUT_DIR = "./results/plots"

METRICS = ["TPR", "FPR", "Precision", "F1"]
MODELS = ["simpleModel", "calibratedModel", "complexModel"]
APPROACHES = ["detection", "accuracy", "mapping"]
CELL_TYPES = ["Normal", "Tumour"]

# Refined color palettes
COLOR_MODEL = {
    "simpleModel":      "#D63031",       # Bold red
    "calibratedModel":  "#0984E3",       # Bold blue
    "complexModel":   "#00B894",       # Bold green
}

COLOR_APPROACH = {
    "detection":  "#E17055",             # Coral
    "accuracy":   "#6C5CE7",             # Violet
    "mapping":    "#00CEC9",             # Cyan
}

COLOR_CELL = {
    "Normal": "#A29BFE",                 # Soft purple
    "Tumour": "#FDCB6E",                 # Warm gold
}

# ── Data loading ────────────────────────────────────────────────────────────────
def load_csv(path):
    """Load per_replicate.csv into a list of dicts."""
    rows = []
    try:
        with open(path) as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                rows.append(row)
        print(f"[load] {len(rows)} rows from {path}")
        return rows
    except Exception as e:
        print(f"[error] Cannot load {path}: {e}")
        return []

def convert_numeric(rows):
    """Convert numeric fields to floats."""
    numeric_fields = ["D", "V", "Replicate", "N_Cells"] + [
        f"{m}_{s}" for m in METRICS for s in ["mean", "std"]
    ]
    for row in rows:
        for field in numeric_fields:
            if field in row and row[field]:
                try:
                    row[field] = float(row[field])
                except (ValueError, TypeError):
                    row[field] = None
    return rows

# ── Data aggregation ────────────────────────────────────────────────────────────
def aggregate_by_dv(rows, metric, model=None, approach=None, cell_type=None):
    """Group rows by (D, V) and return all values + aggregates.
    
    Returns:
        {(d, v): {"values": [val1, val2, ...], "mean": float, "std": float, "n": int}}
    """
    filtered = rows
    if model:
        filtered = [r for r in filtered if r.get("Model") == model]
    if approach:
        filtered = [r for r in filtered if r.get("Approach") == approach]
    if cell_type:
        filtered = [r for r in filtered if r.get("CellType") == cell_type]
    
    agg = defaultdict(list)
    for row in filtered:
        d = row.get("D")
        v = row.get("V")
        if not d or not v:
            continue
        key = (int(d), int(v))
        val = row.get(f"{metric}_mean")
        if val is not None:
            agg[key].append(val)
    
    result = {}
    for key, values in agg.items():
        if values:
            result[key] = {
                "values": values,
                "mean": float(np.mean(values)),
                "std": float(np.std(values) if len(values) > 1 else 0),
                "n": len(values),
            }
    return result

# ── Main plotting function ─────────────────────────────────────────────────────
def plot_scatter_enhanced(rows, metric, group_by, group_val, color, output_dir):
    """Create a single scatter plot with:
      - Light dots for individual replicates (jittered around x)
      - Larger, darker dots for the mean at each D/V
      - Error bars showing ±1 std
    
    Args:
        group_by: "model" | "approach" | "celltype"
        group_val: value of that group (e.g., "simpleModel")
        color: color hex for this group
    """
    # Filter data
    if group_by == "model":
        filtered = [r for r in rows if r.get("Model") == group_val]
    elif group_by == "approach":
        filtered = [r for r in rows if r.get("Approach") == group_val]
    else:  # celltype
        filtered = [r for r in rows if r.get("CellType") == group_val]
    
    # Aggregate by D/V
    agg = aggregate_by_dv(filtered, metric)
    if not agg:
        return False
    
    dv_pairs = sorted(agg.keys())
    
    # Create figure
    fig, ax = plt.subplots(figsize=(13, 7))
    
    # Plot individual replicates + means
    for x_idx, dv in enumerate(dv_pairs):
        data = agg[dv]
        values = data["values"]
        mean_val = data["mean"]
        std_val = data["std"]
        n = data["n"]
        
        # Jitter individual replicates (light dots)
        jitter = np.random.normal(0, 0.03, size=len(values))
        x_jitter = x_idx + jitter
        ax.scatter(x_jitter, values, color=color, alpha=0.25, s=60,
                  edgecolors="none", label="Individual replicates" if x_idx == 0 else "")
        
        # Error bar
        ax.errorbar(x_idx, mean_val, yerr=std_val, fmt='none',
                   ecolor=color, elinewidth=2.5, capsize=7, alpha=0.7, zorder=10)
        
        # Mean dot (larger, darker, with black edge)
        ax.scatter(x_idx, mean_val, color=color, alpha=0.95, s=280,
                  edgecolors="black", linewidth=2, zorder=15,
                  label="Mean ± 1 std" if x_idx == 0 else "")
    
    # Styling
    ax.set_xticks(range(len(dv_pairs)))
    ax.set_xticklabels([f"D{d}V{v}" for d, v in dv_pairs], fontsize=11, fontweight="bold")
    ax.set_ylabel(metric, fontsize=13, fontweight="bold")
    ax.set_xlabel("D/V Combination", fontsize=13, fontweight="bold")
    
    title = f"{metric} — {group_by.capitalize()}: {group_val}"
    ax.set_title(title, fontsize=15, fontweight="bold", pad=15)
    
    # Grid and limits
    ax.grid(True, alpha=0.2, axis='y', linestyle='-', linewidth=0.5)
    ax.set_axisbelow(True)
    ax.set_ylim([-0.05, 1.1])
    
    # Legend
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles[:2], labels[:2], loc="upper left", fontsize=10,
                 framealpha=0.95, edgecolor="black", fancybox=True)
    
    # Spines
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    
    ax.tick_params(axis='both', labelsize=10)
    
    plt.tight_layout()
    
    # Save
    sanitized_val = group_val.replace("/", "_")
    filename = f"scatter_{metric}_{group_by}_{sanitized_val}.png"
    path = os.path.join(output_dir, filename)
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    
    return True

# ── Comparison panels ──────────────────────────────────────────────────────────
def plot_comparison_panel(rows, metric, output_dir):
    """2×2 panel showing all Models, Approaches, CellTypes, and Approach×CellType."""
    fig, axes = plt.subplots(2, 2, figsize=(18, 13))
    fig.suptitle(f"Metric: {metric.upper()}", fontsize=18, fontweight="bold", y=0.995)
    
    # ── Panel 1: By Model ──────────────────────────────────────────────────────
    ax = axes[0, 0]
    dv_pairs_ref = None
    for model in MODELS:
        agg = aggregate_by_dv(rows, metric, model=model)
        if not agg:
            continue
        dv_pairs = sorted(agg.keys())
        if dv_pairs_ref is None:
            dv_pairs_ref = dv_pairs
        x_pos = np.arange(len(dv_pairs))
        means = [agg[dv]["mean"] for dv in dv_pairs]
        stds = [agg[dv]["std"] for dv in dv_pairs]
        
        ax.errorbar(x_pos, means, yerr=stds, fmt='o', markersize=9,
                   label=model, color=COLOR_MODEL[model], alpha=0.8,
                   capsize=6, linewidth=2.5, elinewidth=2)
    
    ax.set_xticks(x_pos)
    ax.set_xticklabels([f"D{d}V{v}" for d, v in dv_pairs], fontsize=10, fontweight="bold")
    ax.set_ylabel(metric, fontsize=11, fontweight="bold")
    ax.set_title("By Model", fontsize=13, fontweight="bold")
    ax.legend(loc="best", fontsize=10, framealpha=0.95)
    ax.grid(True, alpha=0.25, axis='y')
    ax.set_axisbelow(True)
    ax.set_ylim([-0.05, 1.1])
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    
    # ── Panel 2: By Approach ───────────────────────────────────────────────────
    ax = axes[0, 1]
    for approach in APPROACHES:
        agg = aggregate_by_dv(rows, metric, approach=approach)
        if not agg:
            continue
        dv_pairs = sorted(agg.keys())
        x_pos = np.arange(len(dv_pairs))
        means = [agg[dv]["mean"] for dv in dv_pairs]
        stds = [agg[dv]["std"] for dv in dv_pairs]
        
        ax.errorbar(x_pos, means, yerr=stds, fmt='s', markersize=9,
                   label=approach, color=COLOR_APPROACH[approach], alpha=0.8,
                   capsize=6, linewidth=2.5, elinewidth=2)
    
    ax.set_xticks(x_pos)
    ax.set_xticklabels([f"D{d}V{v}" for d, v in dv_pairs], fontsize=10, fontweight="bold")
    ax.set_ylabel(metric, fontsize=11, fontweight="bold")
    ax.set_title("By Approach", fontsize=13, fontweight="bold")
    ax.legend(loc="best", fontsize=10, framealpha=0.95)
    ax.grid(True, alpha=0.25, axis='y')
    ax.set_axisbelow(True)
    ax.set_ylim([-0.05, 1.1])
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    
    # ── Panel 3: By CellType ───────────────────────────────────────────────────
    ax = axes[1, 0]
    for cell_type in CELL_TYPES:
        agg = aggregate_by_dv(rows, metric, cell_type=cell_type)
        if not agg:
            continue
        dv_pairs = sorted(agg.keys())
        x_pos = np.arange(len(dv_pairs))
        means = [agg[dv]["mean"] for dv in dv_pairs]
        stds = [agg[dv]["std"] for dv in dv_pairs]
        
        ax.errorbar(x_pos, means, yerr=stds, fmt='^', markersize=9,
                   label=cell_type, color=COLOR_CELL[cell_type], alpha=0.8,
                   capsize=6, linewidth=2.5, elinewidth=2)
    
    ax.set_xticks(x_pos)
    ax.set_xticklabels([f"D{d}V{v}" for d, v in dv_pairs], fontsize=10, fontweight="bold")
    ax.set_ylabel(metric, fontsize=11, fontweight="bold")
    ax.set_title("By CellType", fontsize=13, fontweight="bold")
    ax.legend(loc="best", fontsize=10, framealpha=0.95)
    ax.grid(True, alpha=0.25, axis='y')
    ax.set_axisbelow(True)
    ax.set_ylim([-0.05, 1.1])
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    
    # ── Panel 4: By Approach × CellType (grouped) ──────────────────────────────
    ax = axes[1, 1]
    width = 0.12
    offset_base = np.arange(len(dv_pairs))
    offsets = [-1.5, -0.5, 0.5, 1.5, 2.5, 3.5]
    
    for idx, (approach, cell_type) in enumerate([
        (a, c) for a in APPROACHES for c in CELL_TYPES
    ]):
        agg = aggregate_by_dv(rows, metric, approach=approach, cell_type=cell_type)
        if not agg:
            continue
        dv_pairs = sorted(agg.keys())
        means = [agg[dv]["mean"] for dv in dv_pairs]
        stds = [agg[dv]["std"] for dv in dv_pairs]
        
        x_pos = offset_base + (offsets[idx % len(offsets)] * width / 1.5)
        marker = ['o', 's', '^'][APPROACHES.index(approach)]
        label = f"{approach[:3]}/{cell_type[:3]}"
        
        ax.errorbar(x_pos, means, yerr=stds, fmt=marker, markersize=6,
                   label=label, alpha=0.65, capsize=3, linewidth=1.5, elinewidth=1)
    
    ax.set_xticks(offset_base)
    ax.set_xticklabels([f"D{d}V{v}" for d, v in dv_pairs], fontsize=10, fontweight="bold")
    ax.set_ylabel(metric, fontsize=11, fontweight="bold")
    ax.set_title("By Approach × CellType", fontsize=13, fontweight="bold")
    ax.legend(loc="best", fontsize=9, ncol=2, framealpha=0.95)
    ax.grid(True, alpha=0.25, axis='y')
    ax.set_axisbelow(True)
    ax.set_ylim([-0.05, 1.1])
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    
    plt.tight_layout()
    path = os.path.join(output_dir, f"panel_{metric}.png")
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"[plot] {path}")

# ── Entry point ────────────────────────────────────────────────────────────────
def main():
    if not os.path.isfile(CSV_PATH):
        print(f"[error] {CSV_PATH} not found.\n  → Run analyseData.py first to generate results/per_replicate.csv")
        sys.exit(1)
    
    print(f"\n{'='*70}")
    print(f"  Enhanced Scatter Plotting from: {CSV_PATH}")
    print(f"{'='*70}\n")
    
    rows = load_csv(CSV_PATH)
    rows = convert_numeric(rows)
    
    if not rows:
        print("[error] No data loaded.")
        sys.exit(1)
    
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Generate comparison panels (one per metric)
    print(f"\n[1/2] Generating comparison panels...")
    for metric in METRICS:
        plot_comparison_panel(rows, metric, OUTPUT_DIR)
    
    # Generate individual scatter plots (one per metric × grouping)
    print(f"\n[2/2] Generating individual scatter plots...")
    total = 0
    
    # By Model
    for metric in METRICS:
        for model in MODELS:
            if plot_scatter_enhanced(rows, metric, "model", model, COLOR_MODEL[model], OUTPUT_DIR):
                total += 1
    
    # By Approach
    for metric in METRICS:
        for approach in APPROACHES:
            if plot_scatter_enhanced(rows, metric, "approach", approach, COLOR_APPROACH[approach], OUTPUT_DIR):
                total += 1
    
    # By CellType
    for metric in METRICS:
        for cell_type in CELL_TYPES:
            if plot_scatter_enhanced(rows, metric, "celltype", cell_type, COLOR_CELL[cell_type], OUTPUT_DIR):
                total += 1
    
    print(f"\n{'='*70}")
    print(f"[done] {len(METRICS)} panels + {total} scatter plots → {OUTPUT_DIR}")
    print(f"{'='*70}\n")

if __name__ == "__main__":
    main()
