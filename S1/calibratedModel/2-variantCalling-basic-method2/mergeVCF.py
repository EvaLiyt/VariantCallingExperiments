#!/usr/bin/env python3
import os
import re
import sys
from typing import Dict, List, Tuple, Optional

BASE_DIR = "../2-variantCalling"
NO_DIR   = os.path.join(BASE_DIR, "D1V1", "merge")  # no-filter baseline

MISSING_FULL = "./.:.:.:.:."

FNAME_RE = re.compile(r"^variants_r(\d+)_d(\d+)_v(\d+)\.vcf$")

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

def corroborated_alleles(row: List[str], alts_with: List[str], sample_cols: List[int]) -> set:
    """Alleles actually variant-called by a strongly-evidenced cell at this site.

    An allele is corroborated only if it appears at index >= 1 in some
    non-missing (i.e. strongly called) cell's genotype in the strict row. An
    allele merely listed in the ALT column but only ever called 0/0 is NOT
    corroborated, so weak baseline calls for it must not be rescued.
    """
    called_idx = set()
    for ci in sample_cols:
        if ci >= len(row):
            continue
        gt = get_gt(row[ci])
        if gt is None:
            continue
        try:
            _, idxs = split_gt(gt)
        except ValueError:
            continue
        called_idx.update(i for i in idxs if i >= 1)
    return {alts_with[i - 1] for i in called_idx if i - 1 < len(alts_with)}

def rescue_corroborated(gt_no: str, alts_no: List[str], alts_with: List[str],
                        corroborated: set) -> str:
    """Map a lenient (baseline) GT onto the strict row's alleles.

    ``alts_with`` is the strict row's ALT list (the output index space, left
    unchanged). ``corroborated`` is the subset of those alleles that some
    strongly-called cell variant-called at this site. REF stays REF, a
    corroborated baseline alt maps to its index in the strict row, and any other
    baseline alt (uncorroborated, or absent from the strict row) collapses to
    REF. No allele is ever added to the row.
    """
    sep, idxs = split_gt(gt_no)
    new_idxs: List[int] = []
    for i in idxs:
        if i <= 0 or i - 1 >= len(alts_no):
            new_idxs.append(0)
            continue
        allele = alts_no[i - 1]
        new_idxs.append(alts_with.index(allele) + 1 if allele in corroborated else 0)
    return join_gt(sep, new_idxs)

def gt_allele_key(gt: str, ref: str, alts: List[str]):
    """Resolve a GT's indices to actual bases for order-independent comparison.

    The two files order their ALT columns independently, so comparing raw index
    strings is unreliable; this maps each index to its base. Falls back to the
    raw GT string if it cannot be parsed (e.g. partial missingness like ``0/.``),
    so it never raises on odd genotypes.
    """
    try:
        _, idxs = split_gt(gt)
    except ValueError:
        return ("raw", gt)
    bases: List[str] = []
    for i in idxs:
        if i == 0:
            bases.append(ref)
        elif 0 < i <= len(alts):
            bases.append(alts[i - 1])
        else:
            bases.append(f"?{i}")
    return ("alleles", tuple(sorted(bases)))

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

# ── per-folder processing ──────────────────────────────────────────────────────

def process_folder(with_dir: str, no_dir: str, out_dir: str) -> None:
    os.makedirs(out_dir, exist_ok=True)

    with_files = sorted([f for f in os.listdir(with_dir) if f.endswith(".vcf")])
    if not with_files:
        print(f"  [WARN] no .vcf files in {with_dir}", file=sys.stderr)
        return

    for fname in with_files:
        m = FNAME_RE.match(fname)
        if not m:
            print(f"  [WARN] {fname} does not match expected pattern, skipping", file=sys.stderr)
            continue

        r_val = m.group(1)
        # always point to D1V1's version of this r value
        no_fname = f"variants_r{r_val}_d1_v1.vcf"

        with_path = os.path.join(with_dir, fname)
        no_path   = os.path.join(no_dir,   no_fname)
        out_path  = os.path.join(out_dir,  fname)
        diff_path = os.path.join(out_dir,  fname.replace(".vcf", ".diff.txt"))

        if not os.path.exists(no_path):
            print(f"  [WARN] baseline {no_fname} not found in {no_dir}, skipping", file=sys.stderr)
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

                # the with-filter row's ALTs are the strong-evidence (corroborated)
                # alleles at this site; the baseline orders its ALTs independently
                alts_with = parse_alts(w[4])
                alts_no   = parse_alts(n[4])
                ref_n     = n[3] if len(n) > 3 else ref_w

                # diff: log GT mismatches where both sides are non-missing,
                # comparing resolved alleles rather than raw index strings
                for sname in with_samples:
                    wi = with_name_to_idx[sname]
                    if sname not in no_name_to_idx:
                        continue
                    ni = no_name_to_idx[sname]
                    if ni >= len(n):
                        continue
                    gt_w = get_gt(w[wi])
                    gt_n = get_gt(n[ni])
                    if gt_w is None or gt_n is None:
                        continue
                    if gt_allele_key(gt_w, ref_w, alts_with) != gt_allele_key(gt_n, ref_n, alts_no):
                        diff_records.append(
                            f"{fname}\t{chrom}\t{pos}\t{sname}\t{gt_w}\t{gt_n}"
                        )

                # fill: rescue missing-in-with cells from the baseline, but only
                # for alleles a strongly-called cell actually variant-called here
                # (not merely listed in ALT), and never adding alleles to the row
                corroborated = corroborated_alleles(
                    w, alts_with, [with_name_to_idx[s] for s in with_samples])

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
                        new_gt = rescue_corroborated(gt_no, alts_no, alts_with, corroborated)
                        if ":" in cell_n:
                            w[wi] = new_gt + ":" + cell_n.split(":", 1)[1]
                        else:
                            w[wi] = new_gt

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

# ── main ───────────────────────────────────────────────────────────────────────

def main():
    dvv_pattern = re.compile(r"^D(\d+)V(\d+)$")

    all_folders = sorted([
        d for d in os.listdir(BASE_DIR)
        if dvv_pattern.match(d)
        and d != "D1V1"
        and os.path.isdir(os.path.join(BASE_DIR, d))
    ])

    if not all_folders:
        print("No D*V* folders found (excluding D1V1). Nothing to do.", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(all_folders)} folder(s): {', '.join(all_folders)}")
    print(f"No-filter baseline : {NO_DIR}\n")

    for folder in all_folders:
        with_dir = os.path.join(BASE_DIR, folder, "merge")
        out_dir  = os.path.join(folder, "merge")          # ./D*V*/merge in CWD

        if not os.path.isdir(with_dir):
            print(f"[WARN] {folder}/merge/ does not exist, skipping", file=sys.stderr)
            continue

        print(f"Processing {BASE_DIR}/{folder}/merge/ -> ./{folder}/merge/")
        process_folder(with_dir, NO_DIR, out_dir)
        print()

if __name__ == "__main__":
    main()