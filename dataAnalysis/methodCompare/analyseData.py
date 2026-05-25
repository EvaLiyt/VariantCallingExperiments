#!/usr/bin/env python3
"""
analyseData.py (rewritten)
──────────────────────────
Run from: fakeFolder/dataAnalysis/methodCompare/

Dynamically reads D (min_depth) and V (min_var) from:
  BASE/{model}/1-simulation/pileup/r{r}/parameters_r{r}.txt

Instead of hard-coded D/V combinations, discovers the actual parameters from
each replicate. This means ~100 replicates total, ~6-7 per D/V combination.

Three approaches across simpleModel, realisticModel, calibratedModel:

DETECTION  (VCF binary, M1/M2/M4)
──────────────────────────────────
Was the variant detected at all?
  TP = variant site, VCF non-ref call
  FN = variant site, VCF missing / ./. / 0/0
  FP = non-variant site, VCF non-ref call
  TN = non-variant site, VCF missing / ./. / 0/0

ACCURACY  (VCF genotype only, M1/M2/M4)
────────────────────────────────────────
Was the variant-caller genotype correct?  No alignment fallback.
  Variant site:
    VCF non-ref, called == true_gt → TP
    VCF non-ref, called != true_gt → FP + FN  (GA4GH: wrong call hurts both)
    VCF 0/0 or missing/./. → FN
  Non-variant site:
    VCF non-ref → FP
    VCF 0/0 or missing/./. → TN

MAPPING  (FASTA alignment only, M1/M2/M3/M4)
─────────────────────────────────────────────
Was the mapper's called alignment genotype correct?  Independent of VCF.
  Variant site:
    FASTA call == true_gt → TP
    FASTA call != true_gt → FP + FN  (GA4GH: wrong call hurts both)
    FASTA missing / ambiguous → FN
  Non-variant site:
    FASTA call == root_gt or missing → TN
    FASTA call != root_gt → FP

OUTPUT
──────
  results/per_cell.csv        — one row per (model, method, D, V, rep, cell, approach)
  results/per_replicate.csv   — avg TPR/FPR/Precision/F1 per (model, method, D, V, rep,
                                 approach, cell_type∈{Normal,Tumour})
"""

import os
import re
import csv
import math
from collections import defaultdict

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE       = "../.."  
OUTPUT_DIR = "./results"

MODELS   = ["simpleModel", "calibratedModel", "realisticModel"]

# ── Genotype helpers ───────────────────────────────────────────────────────────
# Normalise the 16 diploid-hex states to the 10 unordered IUPAC equivalents.
# Symmetric pairs collapse: '1'/'4' → 'm', '2'/'8' → 'r', etc.
GENOTYPE_MAP = {
    '1': 'm', '4': 'm',   # AC / CA
    '2': 'r', '8': 'r',   # AG / GA
    '3': 'w', 'c': 'w',   # AT / TA
    '6': 's', '9': 's',   # CG / GC
    '7': 'y', 'd': 'y',   # CT / TC
    'b': 'k', 'e': 'k',   # GT / TG
}

REF_SKIP = {'-', 'n', '?', '.'}

def norm(g):
    """Canonical form for one genotype character."""
    g = str(g).lower()
    return GENOTYPE_MAP.get(g, g)

def geno_match(a, b):
    return norm(a) == norm(b)

def get_genotype(a1, a2):
    """Convert two nucleotide bases to a diploid-hex char."""
    pair_map = {
        ("a", "a"): "0", ("a", "c"): "1", ("a", "g"): "2", ("a", "t"): "3",
        ("c", "a"): "4", ("c", "c"): "5", ("c", "g"): "6", ("c", "t"): "7",
        ("g", "a"): "8", ("g", "c"): "9", ("g", "g"): "a", ("g", "t"): "b",
        ("t", "a"): "c", ("t", "c"): "d", ("t", "g"): "e", ("t", "t"): "f",
    }
    return pair_map.get((a1.lower(), a2.lower()))

# ── Cell classification ────────────────────────────────────────────────────────
def is_tumour(taxon, model):
    if model == "simpleModel":
        return True   # all cells are tumour in the simple model
    return "clade_" in taxon.lower()

# ── Taxon normalisation ────────────────────────────────────────────────────────
_STRIP_SUFFIX = re.compile(r'_r\d+.*$', re.IGNORECASE)

def norm_taxon(name):
    return _STRIP_SUFFIX.sub('', name).lower().strip()

# ── NEXUS parser ───────────────────────────────────────────────────────────────
def parse_nexus(path):
    """Return {taxon: sequence_string} from a nucleotideDiploid16 nexus file."""
    sequences  = {}
    order      = []
    in_matrix  = False
    interleaved = False
    try:
        with open(path) as fh:
            raw = fh.read()
    except Exception as e:
        print(f"  [nexus] Cannot read {path}: {e}")
        return sequences

    if re.search(r'interleave\s*=\s*yes', raw, re.IGNORECASE):
        interleaved = True

    for line in raw.splitlines():
        stripped = re.sub(r'\[.*?\]', '', line.strip()).strip()
        if not stripped:
            continue
        if re.match(r'^matrix\b', stripped, re.IGNORECASE):
            in_matrix = True
            rest = re.sub(r'^matrix\s*', '', stripped, flags=re.IGNORECASE).strip()
            if rest and rest != ';':
                _nexus_add(rest, sequences, order)
            continue
        if not in_matrix:
            continue
        if stripped == ';':
            if not interleaved:
                break
            in_matrix = False
            continue
        _nexus_add(stripped, sequences, order)
    return sequences

def _nexus_add(stripped, sequences, order):
    parts = stripped.split()
    if len(parts) < 2:
        return
    taxon = parts[0]
    seq   = ''.join(parts[1:]).rstrip(';')   # strip trailing ';' from inline matrix rows
    if taxon not in sequences:
        sequences[taxon] = seq
        order.append(taxon)
    else:
        sequences[taxon] += seq

# ── FASTA parser ───────────────────────────────────────────────────────────────
def parse_fasta(path):
    """Return {norm_taxon: sequence_string} from a FASTA alignment file."""
    seqs         = {}
    current_name = None
    current_seq  = []
    try:
        with open(path) as fh:
            for line in fh:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    if current_name is not None:
                        seqs[current_name] = ''.join(current_seq).lower()
                    current_name = norm_taxon(line[1:].split()[0])
                    current_seq  = []
                else:
                    current_seq.append(line.strip())
        if current_name is not None:
            seqs[current_name] = ''.join(current_seq).lower()
    except Exception as e:
        print(f"  [fasta] Cannot read {path}: {e}")
    return seqs

# ── Root & true-alignment loaders ──────────────────────────────────────────────
def get_root_sequence(model):
    root_dir = os.path.join(BASE, model, "1-simulation", "root")
    if not os.path.isdir(root_dir):
        print(f"  [root] Not found: {root_dir}")
        return None
    for fn in sorted(os.listdir(root_dir)):
        if fn.endswith(".nexus") and not fn.startswith('.'):
            seqs = parse_nexus(os.path.join(root_dir, fn))
            if seqs:
                seq = list(seqs.values())[0]
                print(f"  [root] {model}: {len(seq)} bp from {fn}")
                return seq
    print(f"  [root] No nexus in {root_dir}")
    return None

def index_true_alignments(model):
    nex_dir = os.path.join(BASE, model, "1-simulation", "nexus")
    r_map   = {}
    if not os.path.isdir(nex_dir):
        print(f"  [true aln] Not found: {nex_dir}")
        return r_map
    for fn in sorted(os.listdir(nex_dir)):
        if fn.startswith('.') or not fn.endswith('.nexus'):
            continue
        m = re.search(r'r(\d+)', fn)
        if m:
            r_map[int(m.group(1))] = os.path.join(nex_dir, fn)
    print(f"  [true aln] {model}: {len(r_map)} replicate nexus files")
    return r_map

# ── Parameter file reader ──────────────────────────────────────────────────────
def read_parameters(model, r):
    """Read min_depth and min_var from parameters_r{r}.txt.
    
    Returns: (D, V) or (None, None) if not found.
    """
    param_path = os.path.join(BASE, model, "1-simulation", "pileup", f"r{r}", f"parameters_r{r}.txt")
    if not os.path.isfile(param_path):
        return None, None
    
    try:
        with open(param_path) as fh:
            content = fh.read()
        
        # Parse min_depth: x and min_var: x
        d_match = re.search(r'min_depth\s*:\s*(\d+)', content)
        v_match = re.search(r'min_var\s*:\s*(\d+)', content)
        
        if d_match and v_match:
            d = int(d_match.group(1))
            v = int(v_match.group(1))
            return d, v
    except Exception as e:
        print(f"  [param] Error reading {param_path}: {e}")
    
    return None, None

# ── Method configuration ───────────────────────────────────────────────────────
METHODS = {
    "basic_M1": {
        "type": "basic_dv",
        "dir": "2-variantCalling",
        "do_det": True,
        "do_acc": True,
        "do_map": True
    },
}

# ── Path helpers ───────────────────────────────────────────────────────────────
def _method_root(model, cfg):
    """Construct the root directory for a method's output."""
    return os.path.join(BASE, model, cfg["dir"])

def _vcf_name(r):
    """Return the expected VCF filename."""
    return f"variants_r{r}.vcf"

def _fasta_name(r):
    return f"alignment_r{r}.fasta"

def discover_replicates(model, cfg, true_aln_index):
    """Discover all replicates with both VCF and true alignment data.
    
    For each replicate, read its D/V from parameters.txt.
    Returns:
        {r: {"d": d, "v": v, "vcf": vcf_path, "aln": aln_path, "true_aln": true_aln_path}}
    """
    mroot = os.path.join(BASE, model, cfg["dir"])
    merge_dir = os.path.join(mroot, "merge")
    aln_dir = os.path.join(mroot, "alignment")
    
    if not os.path.isdir(merge_dir):
        print(f"  [discover] merge dir not found: {merge_dir}")
        return {}
    
    found = {}
    
    # Extract replicate numbers from VCF filenames
    for fn in os.listdir(merge_dir):
        if not fn.endswith('.vcf'):
            continue
        
        m = re.search(r'r(\d+)', fn)
        if not m:
            continue
        
        r = int(m.group(1))
        
        # Only include if we have true alignment for this replicate
        if r not in true_aln_index:
            continue
        
        # Read D and V from parameters file
        d, v = read_parameters(model, r)
        if d is None or v is None:
            continue
        
        vcf_path = os.path.join(merge_dir, fn)
        aln_path = os.path.join(aln_dir, _fasta_name(r)) if os.path.isdir(aln_dir) else None
        
        found[r] = {
            "d": d,
            "v": v,
            "vcf": vcf_path,
            "aln": aln_path,
            "true_aln": true_aln_index[r],
        }
    
    return found

# ── VCF parser ─────────────────────────────────────────────────────────────────
def parse_vcf(vcf_path):
    """Parse a merged VCF.

    Returns:
        positive_sites  dict[cell][pos_str] = genotype_char   (non-ref calls)
        negative_sites  dict[cell][pos_str] = genotype_char   (0/0 calls)
    """
    positive_sites = defaultdict(dict)
    negative_sites = defaultdict(dict)
    sample_names   = []
    try:
        with open(vcf_path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith("#CHROM"):
                    sample_names = [norm_taxon(c) for c in line.split("\t")[9:]]
                    continue
                if line.startswith("#") or not line.strip():
                    continue
                cols = line.split("\t")
                if len(cols) < 9:
                    continue
                pos        = cols[1]
                ref        = cols[3].lower()
                alt        = cols[4].lower().split(",")
                fmt_fields = cols[8].split(":")
                gt_idx     = fmt_fields.index("GT") if "GT" in fmt_fields else 0
                for i, cell in enumerate(sample_names):
                    raw = cols[9 + i] if (9 + i) < len(cols) else "."
                    gt  = raw.split(":")[gt_idx] if raw not in (".", "") else "./."
                    gt  = gt.replace("|", "/")
                    if gt in ("./.", "."):
                        continue
                    parts = gt.split("/")
                    if len(parts) != 2 or "." in parts:
                        continue
                    a1, a2   = parts
                    base1    = ref if a1 == "0" else alt[int(a1) - 1] if int(a1) - 1 < len(alt) else ref
                    base2    = ref if a2 == "0" else alt[int(a2) - 1] if int(a2) - 1 < len(alt) else ref
                    genotype = get_genotype(base1, base2)
                    if genotype is None:
                        continue
                    if a1 == "0" and a2 == "0":
                        negative_sites[cell][pos] = genotype
                    else:
                        positive_sites[cell][pos] = genotype
    except Exception as e:
        print(f"  [vcf] Error: {vcf_path}: {e}")
    return dict(positive_sites), dict(negative_sites)

# ── Core analysis ──────────────────────────────────────────────────────────────

def _variant_sets(root_seq, true_seq):
    """Return (true_variants, true_non_variants) as {pos_str: genotype_char} dicts."""
    max_len = min(len(root_seq), len(true_seq))
    var    = {}
    nonvar = {}
    for i in range(max_len):
        rg = root_seq[i]
        tg = true_seq[i]
        if norm(rg) in REF_SKIP or norm(tg) in REF_SKIP:
            continue
        pos = str(i + 1)
        if geno_match(rg, tg):
            nonvar[pos] = tg
        else:
            var[pos] = tg
    return var, nonvar


def analyse_detection(true_seqs, root_seq, positive_sites, negative_sites):
    """Detection: VCF binary, no alignment fallback.  M1/M2/M4.

    TP = variant site, VCF non-ref call
    FN = variant site, VCF missing / ./. / 0/0
    FP = non-variant site, VCF non-ref call
    TN = non-variant site, VCF missing / ./. / 0/0
    """
    per_taxon = {}
    for taxon, true_seq in true_seqs.items():
        nk       = norm_taxon(taxon)
        cell_pos = positive_sites.get(nk, {})

        var, nonvar = _variant_sets(root_seq, true_seq)

        tp = sum(1 for pos in var    if pos in cell_pos)
        fn = sum(1 for pos in var    if pos not in cell_pos)
        fp = sum(1 for pos in cell_pos if pos not in var)
        tn = sum(1 for pos in nonvar if pos not in cell_pos)

        per_taxon[taxon] = {"TP": tp, "FP": fp, "FN": fn, "TN": tn}
    return per_taxon


def analyse_accuracy(true_seqs, root_seq, positive_sites, negative_sites):
    """Accuracy: VCF genotype only, no alignment fallback.  M1/M2/M4.

    Variant site:
      VCF non-ref, called == true_gt  → TP
      VCF non-ref, called != true_gt  → FP + FN  (GA4GH)
      VCF 0/0 or missing/./.          → FN
    Non-variant site:
      VCF non-ref                     → FP
      VCF 0/0 or missing/./.          → TN
    """
    per_taxon = {}
    unmatched = []
    for taxon, true_seq in true_seqs.items():
        nk       = norm_taxon(taxon)
        cell_pos = positive_sites.get(nk, {})
        cell_neg = negative_sites.get(nk, {})

        if nk not in positive_sites and nk not in negative_sites:
            unmatched.append(taxon)

        var, nonvar = _variant_sets(root_seq, true_seq)

        tp = fp = fn = tn = 0

        for pos, true_g in var.items():
            if pos in cell_pos:
                if geno_match(cell_pos[pos], true_g):
                    tp += 1
                else:
                    fp += 1
                    fn += 1   # GA4GH: wrong genotype hurts both precision and recall
            else:
                fn += 1   # 0/0 or not in VCF → missed variant

        for pos in nonvar:
            if pos in cell_pos:
                fp += 1   # non-ref call at non-variant site
            else:
                tn += 1   # 0/0 or not in VCF → correctly not called

        per_taxon[taxon] = {"TP": tp, "FP": fp, "FN": fn, "TN": tn}

    if unmatched:
        print(f"  [warn] {len(unmatched)} taxa not in VCF: "
              f"{unmatched[:5]}" + (" ..." if len(unmatched) > 5 else ""))
    return per_taxon


def analyse_mapping(true_seqs, root_seq, called_seqs):
    """Mapping: compare called alignment (FASTA) to true alignment (nexus).  M1/M2/M3/M4.

    The comparison is purely called vs true — root is not used in scoring,
    only in _variant_sets to classify each position as variant or non-variant.

    Variant site   (called[i] vs true[i]):
      FASTA == true → TP
      FASTA != true → FP + FN  (GA4GH: wrong call hurts both)
      FASTA missing/? → FN

    Non-variant site (true[i] == root[i] by definition):
      FASTA == true → TN
      FASTA != true → FP
      FASTA missing/? → TN
    """
    per_taxon = {}
    for taxon, true_seq in true_seqs.items():
        nk      = norm_taxon(taxon)
        aln_seq = called_seqs.get(nk, "")

        var, nonvar = _variant_sets(root_seq, true_seq)

        tp = fp = fn = tn = 0

        for pos, true_g in var.items():
            i = int(pos) - 1
            if not aln_seq or i >= len(aln_seq) or norm(aln_seq[i]) in REF_SKIP:
                fn += 1
                continue
            if geno_match(aln_seq[i], true_g):  tp += 1
            else:
                fp += 1   # GA4GH: wrong call hurts both precision and recall
                fn += 1

        for pos, true_g in nonvar.items():
            i = int(pos) - 1
            if not aln_seq or i >= len(aln_seq) or norm(aln_seq[i]) in REF_SKIP:
                tn += 1
                continue
            if geno_match(aln_seq[i], true_g):  tn += 1
            else:                                fp += 1

        per_taxon[taxon] = {"TP": tp, "FP": fp, "FN": fn, "TN": tn}
    return per_taxon

# ── Metrics ────────────────────────────────────────────────────────────────────
METRIC_KEYS = ["TPR", "FPR", "Precision", "F1"]

def calc_metrics(counts):
    tp = counts.get("TP", 0)
    fp = counts.get("FP", 0)
    fn = counts.get("FN", 0)
    tn = counts.get("TN", 0)
    tpr  = tp / (tp + fn)          if (tp + fn) > 0        else 0.0
    fpr  = fp / (fp + tn)          if (fp + tn) > 0        else 0.0
    prec = tp / (tp + fp)          if (tp + fp) > 0        else 0.0
    f1   = 2 * tp / (2*tp + fp + fn) if (2*tp + fp + fn) > 0 else 0.0
    return {"TPR": tpr, "FPR": fpr, "Precision": prec, "F1": f1}

def mean_std(values):
    n = len(values)
    if n == 0:
        return 0.0, 0.0
    mu  = sum(values) / n
    std = math.sqrt(sum((x - mu)**2 for x in values) / n) if n > 1 else 0.0
    return mu, std

# ── Per-model processing ───────────────────────────────────────────────────────
def process_model(model, root_seq, true_aln_index):
    cell_rows = []
    rep_rows  = []

    for method_name, cfg in METHODS.items():
        # Discover replicates and their D/V parameters
        rep_dict = discover_replicates(model, cfg, true_aln_index)
        
        if not rep_dict:
            print(f"    [warn] No replicates discovered for {method_name}")
            continue
        
        print(f"    [discover] Found {len(rep_dict)} replicates with D/V parameters")
        
        for r in sorted(rep_dict.keys()):
            rep_info = rep_dict[r]
            d = rep_info["d"]
            v = rep_info["v"]
            vcf_path = rep_info["vcf"]
            aln_path = rep_info["aln"]
            true_aln_path = rep_info["true_aln"]
            
            true_seqs = parse_nexus(true_aln_path)
            if not true_seqs:
                continue
            
            if not os.path.isfile(vcf_path):
                continue
            
            pos_sites, neg_sites = parse_vcf(vcf_path)
            
            called_seqs = parse_fasta(aln_path) if aln_path and os.path.isfile(aln_path) else {}
            
            # Collect per-approach per-taxon counts
            approach_results = {}
            
            if cfg["do_det"]:
                approach_results["detection"] = analyse_detection(
                    true_seqs, root_seq, pos_sites, neg_sites)
            
            if cfg["do_acc"]:
                approach_results["accuracy"] = analyse_accuracy(
                    true_seqs, root_seq, pos_sites, neg_sites)
            
            if cfg["do_map"]:
                approach_results["mapping"] = analyse_mapping(
                    true_seqs, root_seq, called_seqs)
            
            # Build cell rows and accumulate for per-replicate summary
            agg = defaultdict(lambda: {"Normal": [], "Tumour": []})
            
            for approach, per_taxon in approach_results.items():
                for taxon, counts in per_taxon.items():
                    m         = calc_metrics(counts)
                    cell_type = "Tumour" if is_tumour(taxon, model) else "Normal"
                    cell_rows.append({
                        "Model": model, "Method": method_name,
                        "D": d, "V": v, "Replicate": r,
                        "Taxon": taxon, "CellType": cell_type,
                        "Approach": approach,
                        "TP": counts["TP"], "FP": counts["FP"],
                        "FN": counts["FN"], "TN": counts["TN"],
                        **m,
                    })
                    agg[approach][cell_type].append(m)
            
            # Per-replicate rows
            for approach, by_type in agg.items():
                for cell_type, mlist in by_type.items():
                    if not mlist:
                        continue
                    row = {
                        "Model": model, "Method": method_name,
                        "D": d, "V": v, "Replicate": r,
                        "Approach": approach, "CellType": cell_type,
                        "N_Cells": len(mlist),
                    }
                    for k in METRIC_KEYS:
                        mu, sd = mean_std([m[k] for m in mlist])
                        row[f"{k}_mean"] = mu
                        row[f"{k}_std"]  = sd
                    rep_rows.append(row)
    
    return cell_rows, rep_rows

# ── Output helpers ─────────────────────────────────────────────────────────────
_CELL_FIELDS = (
    ["Model", "Method", "D", "V", "Replicate", "Taxon", "CellType", "Approach",
     "TP", "FP", "FN", "TN"] + METRIC_KEYS
)

def _rep_fields():
    base = ["Model", "Method", "D", "V", "Replicate", "Approach", "CellType", "N_Cells"]
    for k in METRIC_KEYS:
        base += [f"{k}_mean", f"{k}_std"]
    return base

def write_csv(rows, path, fieldnames):
    if not rows:
        print(f"  [output] No rows to write: {path}")
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"  [output] {len(rows)} rows → {path}")

# ── Entry point ────────────────────────────────────────────────────────────────
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    all_cell_rows = []
    all_rep_rows  = []

    for model in MODELS:
        print(f"\n{'='*60}\n  Model: {model}\n{'='*60}")
        root_seq = get_root_sequence(model)
        if root_seq is None:
            continue
        true_aln_index = index_true_alignments(model)
        if not true_aln_index:
            continue
        cell_rows, rep_rows = process_model(model, root_seq, true_aln_index)
        print(f"  → {len(cell_rows)} cell rows, {len(rep_rows)} replicate rows")
        all_cell_rows.extend(cell_rows)
        all_rep_rows.extend(rep_rows)

    write_csv(all_cell_rows, os.path.join(OUTPUT_DIR, "per_cell.csv"), _CELL_FIELDS)
    write_csv(all_rep_rows,  os.path.join(OUTPUT_DIR, "per_replicate.csv"),  _rep_fields())
    print(f"\n[done] {len(all_cell_rows)} cell rows, {len(all_rep_rows)} replicate rows total")


if __name__ == "__main__":
    main()