#!/usr/bin/env python3
"""Tests for 1-mergeVCF.py.

Exercises the real process_folder() rather than a reimplementation. Run with:

    python3 test_mergeVCF.py

Covers the allele-level corroboration rule (a missing-in-strict cell is rescued
from the baseline only for alleles already strongly evidenced at the site, and
the ALT column never gains an uncorroborated allele) and the base-resolved diff
(genotype mismatches compared by allele, not by raw index, so independent ALT
orderings do not create false positives/negatives).
"""
import importlib.util
import os
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("mergevcf", os.path.join(HERE, "mergeVCF.py"))
mv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mv)

_failures = []


def check(name, got, want):
    ok = got == want
    print(f"  [{'OK ' if ok else 'FAIL'}] {name}: got={got!r} want={want!r}")
    if not ok:
        _failures.append(name)


def _run(with_text, no_text, with_name="variants_r0_d4_v4.vcf",
         no_name="variants_r0_d1_v1.vcf"):
    """Run the real process_folder on two in-memory VCFs; return (rows, diffs)."""
    d = tempfile.mkdtemp()
    with_dir = os.path.join(d, "with"); os.makedirs(with_dir)
    no_dir = os.path.join(d, "no"); os.makedirs(no_dir)
    out_dir = os.path.join(d, "out")
    with open(os.path.join(with_dir, with_name), "w") as f:
        f.write(with_text)
    with open(os.path.join(no_dir, no_name), "w") as f:
        f.write(no_text)
    mv.process_folder(with_dir, no_dir, out_dir)
    rows = {}
    with open(os.path.join(out_dir, with_name)) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            rows[p[1]] = p
    diffs = set()
    diff_path = os.path.join(out_dir, with_name.replace(".vcf", ".diff.txt"))
    if os.path.exists(diff_path):
        with open(diff_path) as f:
            for line in f:
                if line.startswith("file"):
                    continue
                c = line.rstrip("\n").split("\t")
                diffs.add((c[2], c[3]))  # (POS, CELL)
    return rows, diffs


HDR = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tc1\tc2\tc3\tc4\n"
)


def test_synthetic():
    print("\n== synthetic (covers corroboration + ALT re-indexing) ==")
    # strict (with-filter): ALT column holds only strongly-evidenced alleles
    with_text = HDR + "\n".join([
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:AD\t0/1:5,5\t./.:.\t./.:.\t0/1:6,4",
        "chr1\t200\t.\tA\tT\t.\t.\t.\tGT:AD\t0/1:5,5\t0/1:5,5\t0/0:9,0\t./.:.",
        "chr1\t300\t.\tA\tT\t.\t.\t.\tGT:AD\t0/0:9,0\t./.:.\t./.:.\t./.:.",
        "chr1\t400\t.\tA\tT\t.\t.\t.\tGT:AD\t0/1:5,5\t./.:.\t./.:.\t./.:.",
        # 500: C is listed in ALT but only ever called 0/0 (a "phantom" allele)
        "chr1\t500\t.\tA\tC\t.\t.\t.\tGT:AD\t0/0:9,0\t./.:.\t./.:.\t./.:.",
    ]) + "\n"
    # baseline (lenient): ALT columns ordered independently, extra artefacts
    no_text = HDR + "\n".join([
        "chr1\t100\t.\tA\tT,G\t.\t.\t.\tGT:AD\t0/1:5,5\t0/2:7,3\t0/1:6,4\t0/1:6,4",
        "chr1\t200\t.\tA\tG,T\t.\t.\t.\tGT:AD\t0/2:4,6\t0/0:9,0\t0/0:9,0\t0/0:9,0",
        "chr1\t300\t.\tA\tT\t.\t.\t.\tGT:AD\t0|0:9,0\t0/0:9,0\t0/0:9,0\t0/0:9,0",
        "chr1\t400\t.\tA\tG,T\t.\t.\t.\tGT:AD\t0/1:5,5\t0/2:4,6\t0/0:9,0\t0/0:9,0",
        "chr1\t500\t.\tA\tC\t.\t.\t.\tGT:AD\t0/0:9,0\t0/1:6,4\t0/0:9,0\t0/0:9,0",
    ]) + "\n"
    rows, diffs = _run(with_text, no_text)

    # ALT columns must stay single-allele (no uncorroborated allele added)
    for pos in ("100", "200", "300", "400", "500"):
        check(f"ALT@{pos} unchanged", rows[pos][4], "C" if pos == "500" else "T")

    # corroboration-aware rescue
    check("100 c1 present untouched", rows["100"][9], "0/1:5,5")
    check("100 c2 baseline G not corroborated -> REF", rows["100"][10], "0/0:7,3")
    check("100 c3 baseline T corroborated -> 0/1", rows["100"][11], "0/1:6,4")
    check("400 c2 baseline idx2=T corroborated -> 0/1", rows["400"][10], "0/1:4,6")
    # phantom allele: C is in ALT but no strong cell variant-called it -> no rescue
    check("500 c2 weak C, allele only ever 0/0 -> REF", rows["500"][10], "0/0:6,4")

    # base-resolved diff
    check("diff count", len(diffs), 2)
    check("real diff (200,c2) logged", ("200", "c2") in diffs, True)
    check("false-negative (400,c1) now caught", ("400", "c1") in diffs, True)
    check("re-index FP (200,c1) not logged", ("200", "c1") in diffs, False)
    check("phasing FP (300,c1) not logged", ("300", "c1") in diffs, False)


def test_real_example():
    """Regression against the committed example VCFs, if present."""
    base = os.path.join(HERE, "..", "2-variantCalling")
    with_dir = os.path.join(base, "D2V1", "merge")
    no_dir = os.path.join(base, "D1V1", "merge")
    name = "variants_r0_d2_v1.vcf"
    if not (os.path.exists(os.path.join(with_dir, name))
            and os.path.exists(os.path.join(no_dir, "variants_r0_d1_v1.vcf"))):
        print("\n== real example: SKIPPED (example VCFs not present) ==")
        return
    print("\n== real example (D2V1 vs D1V1) ==")
    out_dir = tempfile.mkdtemp()
    mv.process_folder(with_dir, no_dir, out_dir)

    in_alts = {}
    with open(os.path.join(with_dir, name)) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            in_alts[p[1]] = set(p[4].split(","))

    out_rows = {}
    violations = 0
    with open(os.path.join(out_dir, name)) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            out_rows[p[1]] = p
            if set(p[4].split(",")) - in_alts.get(p[1], set()):
                violations += 1
    # invariant: no allele is ever introduced beyond the strict input's ALTs
    check("no allele introduced beyond strict input", violations, 0)
    # specific known record: chr1:1081 single-read baseline G must not be minted
    if "1081" in out_rows:
        check("chr1:1081 ALT stays 'C'", out_rows["1081"][4], "C")
    # chr1:124 lists ALT=G but every strong cell is 0/0 (phantom); after fill no
    # cell may be variant-called there
    if "124" in out_rows:
        variant_cells = 0
        for cell in out_rows["124"][9:]:
            gt = mv.get_gt(cell)
            if gt is None:
                continue
            try:
                _, idxs = mv.split_gt(gt)
            except ValueError:
                continue
            if any(i >= 1 for i in idxs):
                variant_cells += 1
        check("chr1:124 phantom G never rescued as a variant", variant_cells, 0)


if __name__ == "__main__":
    test_synthetic()
    test_real_example()
    print()
    if _failures:
        print(f"FAILED: {len(_failures)} check(s): {_failures}")
        raise SystemExit(1)
    print("ALL CHECKS PASSED")
