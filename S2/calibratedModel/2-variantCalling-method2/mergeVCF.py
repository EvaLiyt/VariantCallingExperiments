#!/usr/bin/env python3
import os
import re
import sys
from typing import Dict, List, Tuple, Optional

BASE_DIR = "../2-variantCalling"
NO_DIR = "./merge"  # no-filter baseline

MISSING_FULL = "./.:.:.:.:."

FNAME_RE = re.compile(r"^variants_r(\d+)\.vcf$")

# ── helpers ────────────────────────────────────────────────────────────────────

def is_missing(sample: str) -> bool:
    return (
        sample == MISSING_FULL
        or sample == "."
        or sample.startswith("./.")
        or sample.startswith(".|.")
        or sample.startswith("./")
        or sample.startswith(".|")
    )

def get_gt(sample: str) -> Optional[str]:
    if is_missing(sample):
        return None
    return sample.split(":", 1)[0]

def split_gt(gt: str) -> Tuple[str, List[int]]:
    sep = "/" if "/" in gt else "|"
    return sep, [int(x) for x in gt.split(sep)]

def join_gt(sep: str, idxs: List[int]) -> str:
    return sep.join(str(i) for i in idxs)

def parse_alts(alt: str) -> List[str]:
    return [] if alt in (".", "") else alt.split(",")

def remap_gt_from_no_to_target(gt_no: str, alts_no: List[str], alts_target: List[str]) -> str:
    sep, idxs = split_gt(gt_no)
    new_idxs: List[int] = []
    for i in idxs:
        if i == 0:
            new_idxs.append(0)
            continue
        if i - 1 < 0 or i - 1 >= len(alts_no):
            new_idxs.append(i)
            continue
        allele = alts_no[i - 1]
        if allele not in alts_target:
            alts_target.append(allele)
        new_idxs.append(alts_target.index(allele) + 1)
    return join_gt(sep, new_idxs)

def build_target_alts(alts_with: List[str], alts_no: List[str], gt_no: str) -> List[str]:
    _, idxs = split_gt(gt_no)
    used: List[str] = []
    for i in idxs:
        if i > 0 and (i - 1) < len(alts_no):
            used.append(alts_no[i - 1])
    target = alts_with[:]
    for a in used:
        if a not in target:
            target.append(a)
    return target

def transform_sample_to_target(sample_no: str, alts_no: List[str], alts_target: List[str]) -> str:
    gt_no = get_gt(sample_no)
    if gt_no is None:
        return sample_no
    new_gt = remap_gt_from_no_to_target(gt_no, alts_no, alts_target)
    if ":" in sample_no:
        rest = sample_no.split(":", 1)[1]
        return new_gt + ":" + rest
    return new_gt

def parse_no_vcf_index(no_path: str) -> Tuple[Dict[Tuple[str, str], List[List[str]]], List[str]]:
    index: Dict[Tuple[str, str], List[List[str]]] = {}
    no_samples: List[str] = []
    with open(no_path, "r") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                hdr = line.rstrip("\n").split("\t")
                no_samples = hdr[9:]
                continue
            parts = line.rstrip("\n").split("\t")
            chrom, pos = parts[0], parts[1]
            index.setdefault((chrom, pos), []).append(parts)
    return index, no_samples

def choose_best_no_record(no_records: List[List[str]], ref_w: str) -> Optional[List[str]]:
    if not no_records:
        return None
    for r in no_records:
        if len(r) > 3 and r[3] == ref_w:
            return r
    return no_records[0]

# ── main ───────────────────────────────────────────────────────────────────────

def main():
    """
    Process VCF files in ./merge directory.
    Expects files named: variants_r<N>.vcf
    
    Maps each file to baseline: ./merge/variants_r<N>.vcf
    Outputs merged results to: ./merge/variants_r<N>.vcf (overwrites input)
    """
    with_dir = "./merge"

    if not os.path.isdir(with_dir):
        print(f"[ERROR] {with_dir}/ directory not found", file=sys.stderr)
        sys.exit(1)

    with_files = sorted([f for f in os.listdir(with_dir) if f.endswith(".vcf")])
    if not with_files:
        print(f"[ERROR] no .vcf files in {with_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(with_files)} VCF file(s) in {with_dir}/")
    print(f"No-filter baseline : {NO_DIR}\n")

    for fname in with_files:
        m = FNAME_RE.match(fname)
        if not m:
            print(f"  [WARN] {fname} does not match expected pattern (variants_r<N>.vcf), skipping", file=sys.stderr)
            continue

        r_val = m.group(1)
        # Point to baseline file with same r value
        no_fname = f"variants_r{r_val}.vcf"

        with_path = os.path.join(with_dir, fname)
        no_path   = os.path.join(NO_DIR, no_fname)
        out_fname = f"combined_r{r_val}.vcf"
        out_path  = os.path.join(with_dir, out_fname)
        diff_fname = f"combined_r{r_val}.diff.txt"
        diff_path = os.path.join(with_dir, diff_fname)

        if not os.path.exists(no_path):
            print(f"  [WARN] baseline {no_fname} not found in {NO_DIR}, skipping", file=sys.stderr)
            continue

        no_index, no_samples = parse_no_vcf_index(no_path)
        no_name_to_idx = {name: 9 + i for i, name in enumerate(no_samples)}

        diff_records: List[str] = []

        with open(with_path, "r") as fw, open(out_path, "w") as fo:
            with_samples: List[str] = []
            with_name_to_idx: Dict[str, int] = {}

            for line in fw:
                if line.startswith("##"):
                    fo.write(line)
                    continue

                if line.startswith("#CHROM"):
                    hdr = line.rstrip("\n").split("\t")
                    with_samples = hdr[9:]
                    with_name_to_idx = {name: 9 + i for i, name in enumerate(with_samples)}
                    fo.write(line)
                    continue

                w = line.rstrip("\n").split("\t")
                chrom, pos, ref_w = w[0], w[1], w[3]

                no_records = no_index.get((chrom, pos), [])
                n = choose_best_no_record(no_records, ref_w)
                if n is None:
                    fo.write("\t".join(w) + "\n")
                    continue

                # diff: log GT mismatches where both sides are non-missing
                for sname in with_samples:
                    wi = with_name_to_idx[sname]
                    if sname not in no_name_to_idx:
                        continue
                    ni = no_name_to_idx[sname]
                    if ni >= len(n):
                        continue
                    gt_w = get_gt(w[wi])
                    gt_n = get_gt(n[ni])
                    if gt_w is not None and gt_n is not None and gt_w != gt_n:
                        diff_records.append(
                            f"{fname}\t{chrom}\t{pos}\t{sname}\t{gt_w}\t{gt_n}"
                        )

                # fill: patch missing-in-with from no
                alts_with = parse_alts(w[4])
                alts_no   = parse_alts(n[4])

                for sname in with_samples:
                    wi = with_name_to_idx[sname]
                    if sname not in no_name_to_idx:
                        continue
                    ni = no_name_to_idx[sname]
                    if ni >= len(n):
                        continue
                    cell_w = w[wi]
                    cell_n = n[ni]

                    if is_missing(cell_w) and not is_missing(cell_n):
                        gt_no = get_gt(cell_n)
                        if gt_no is None:
                            continue
                        target_alts = build_target_alts(parse_alts(w[4]), alts_no, gt_no)
                        new_sample  = transform_sample_to_target(cell_n, alts_no, target_alts)
                        w[4]  = ",".join(target_alts) if target_alts else "."
                        w[wi] = new_sample

                fo.write("\t".join(w) + "\n")

        if diff_records:
            with open(diff_path, "w") as fd:
                fd.write("file\tCHROM\tPOS\tCELL\tGT_withFilter\tGT_noFilter\n")
                for r in diff_records:
                    fd.write(r + "\n")
            print(f"  [OK] {fname} <- {no_fname} | merged, diff written ({len(diff_records)} lines)")
        else:
            if os.path.exists(diff_path):
                os.remove(diff_path)
            print(f"  [OK] {fname} <- {no_fname} | merged, no diff")

if __name__ == "__main__":
    main()