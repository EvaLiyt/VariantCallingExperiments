#!/usr/bin/env python3
"""
update_bulk_lphy.py

For each replicate r (0..99), reads parameter values from simulation outputs
and rewrites the corresponding bulk/bulk_r{r}.lphy file with those values.

─── CONFIG ──────────────────────────────────────────────────────────────────
Edit the section below to choose which parameters to trace.
Set a parameter to True to trace it (replace in .lphy), False to skip.
─────────────────────────────────────────────────────────────────────────────
"""

import os
import re
import sys

# ─── USER CONFIG ─────────────────────────────────────────────────────────────

BULK_DIR       = "bulkGenerationScripts"
LPHY_DIR       = "lphy_scripts"        # contains sequenceGeneration_rx_psihealthy.trees
SIM_DATA_DIR   = "simulationData"

# Range of replicates
R_START = 0
R_END   = 99   # inclusive

# Scalar parameters: key must match EXACTLY the column header in the .log file.
# rootAge is special — it is parsed from the .trees file instead of the .log.
TRACE_SCALARS = {
    "lambda":   True,
    "rootRate": True,
    "rootAge":  True,   # special: parsed from lphy_scripts/sequenceGeneration_rx_psihealthy.trees
}

# Vector parameters: each entry defines a prefix and how many indexed columns
# to collect (prefix_0, prefix_1, ..., prefix_{size-1}).
# lphy_var is the variable name to match and replace in the .lphy file.
TRACE_VECTORS = {
    "pi": {
        "active":   True,
        "prefix":   "pi",       # column headers: pi_0, pi_1, ..., pi_15
        "size":     16,
        "lphy_var": "π",        # matches the  π ~  line in the .lphy file
    },
    "rates": {
        "active":   True,
        "prefix":   "rates",    # column headers: rates_0, ..., rates_5
        "size":     6,
        "lphy_var": "rates",
    },
}

# ─── END CONFIG ──────────────────────────────────────────────────────────────


def parse_log(log_path: str, wanted: set) -> dict:
    """
    Read the first two non-comment lines of a BEAST/LPhY .log file.
    Line 1 = tab-separated header.
    Line 2 = tab-separated values (first sample row).

    Returns a dict {column_name: value_string} for columns in `wanted`.
    """
    with open(log_path, "r") as fh:
        lines = [l.rstrip("\n") for l in fh if not l.startswith("#")]

    if len(lines) < 2:
        raise ValueError(f"Log file too short: {log_path}")

    headers = lines[0].split("\t")
    values  = lines[1].split("\t")

    result = {}
    for name in wanted:
        if name in headers:
            result[name] = values[headers.index(name)]
        else:
            raise KeyError(f"Column '{name}' not found in {log_path}")
    return result


def parse_vector_from_log(log_path: str, prefix: str, size: int) -> list[float]:
    """
    Collect prefix_0, prefix_1, ..., prefix_{size-1} from the first data row
    and return them as a list of floats.
    """
    wanted = {f"{prefix}_{i}" for i in range(size)}
    raw = parse_log(log_path, wanted)
    return [float(raw[f"{prefix}_{i}"]) for i in range(size)]


def format_vector_lphy(values: list[float]) -> str:
    """Format a list of floats as a LPhY 1-D array literal: [v0, v1, ...]"""
    return "[" + ", ".join(f"{v:.10g}" for v in values) + "]"


def get_root_age(trees_path: str) -> float:
    """
    Parse a NEXUS .trees file, extract the first Newick tree string,
    and return the root age (= maximum root-to-tip path length).

    Uses [&height=...] node annotations when present (BEAST output),
    otherwise sums branch lengths along the deepest path.
    """
    newick = _extract_first_newick(trees_path)
    return _root_age_from_newick(newick)


def _extract_first_newick(trees_path: str) -> str:
    with open(trees_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if line.lower().startswith("tree "):
                _, _, rest = line.partition("=")
                return re.sub(r"^\s*\[&[^\]]*\]\s*", "", rest.strip())
    raise ValueError(f"No tree found in {trees_path}")


def _root_age_from_newick(newick: str) -> float:
    # Prefer [&height=X] annotation on the root node
    m = re.search(r"\)\s*\[&[^\]]*\bheight\s*=\s*([0-9Ee.+\-]+)", newick)
    if m:
        return float(m.group(1))
    return _max_path_length(newick)


def _max_path_length(newick: str) -> float:
    """Stack-based Newick parser that returns the deepest root-to-tip distance."""
    s = newick.strip().rstrip(";").strip()
    stack: list[float] = [0.0]
    max_depth = 0.0
    i = 0

    while i < len(s):
        c = s[i]

        if c == "(":
            stack.append(stack[-1])
            i += 1

        elif c == ")":
            i += 1
            while i < len(s) and s[i] not in (":", ",", ")", "(", "[", ";"):
                i += 1
            if i < len(s) and s[i] == "[":
                end = s.index("]", i)
                i = end + 1
            bl = 0.0
            if i < len(s) and s[i] == ":":
                i += 1
                m = re.match(r"[0-9Ee.+\-]+", s[i:])
                if m:
                    bl = float(m.group())
                    i += len(m.group())
            child_depth = stack.pop() + bl
            if stack:
                stack[-1] = max(stack[-1], child_depth)
            else:
                max_depth = max(max_depth, child_depth)

        elif c == ",":
            sibling_depth = stack.pop()
            max_depth = max(max_depth, sibling_depth)
            stack.append(stack[-1] if stack else 0.0)
            i += 1

        elif c == "[":
            end = s.index("]", i)
            i = end + 1

        elif c == ":":
            i += 1
            m = re.match(r"[0-9Ee.+\-]+", s[i:])
            if m:
                bl = float(m.group())
                i += len(m.group())
                stack[-1] += bl
            max_depth = max(max_depth, stack[-1])

        else:
            i += 1

    return max(max_depth, max(stack) if stack else 0.0)


def update_lphy(lphy_path: str,
                scalar_replacements: dict,
                vector_replacements: dict) -> None:
    """
    Read bulk_r{x}.lphy, apply all replacements, write it back.

    scalar_replacements: { param_name: value_str }
        e.g. {"lambda": "2.0", "rootAge": "8.34", "rootRate": "0.5"}

    vector_replacements: { lphy_var: formatted_array_str }
        e.g. {"π": "[0.1, 0.2, ...]", "rates": "[0.3, ...]"}

    Rules
    ─────
    • all scalars  → replaces the line  `<param> ~ ...`  or  `<param> = ...`
                     with               `<param> = <value>;`
                     rootAge is only replaced in the model{} block;
                     if it exists in the data{} block it is removed.
    • vector       → replaces the line  `<lphy_var> ~ ...`
                     with               `<lphy_var> = <array>;`
    """
    with open(lphy_path, "r") as fh:
        lines = fh.readlines()

    new_lines = []
    in_data_block = False

    for line in lines:
        stripped = line.rstrip()
        indent = len(line) - len(line.lstrip())

        # ── Track data{} block ────────────────────────────────────────────
        if re.search(r"\bdata\s*\{", stripped):
            in_data_block = True
        if in_data_block and "}" in stripped:
            in_data_block = False

        # ── Drop rootAge if it appears inside the data{} block ────────────
        if in_data_block and re.match(r"\s*rootAge\s*[=~]", stripped):
            continue

        # ── Vector parameters (π, rates, …) ──────────────────────────────
        matched_vec = None
        for lphy_var, arr_str in vector_replacements.items():
            if re.match(rf"\s*{re.escape(lphy_var)}\s*[=~]", stripped):
                matched_vec = (lphy_var, arr_str)
                break

        if matched_vec:
            lphy_var, arr_str = matched_vec
            new_lines.append(" " * indent + f"{lphy_var} = {arr_str};\n")
            continue

        # ── Scalar parameters (lambda, rootRate, rootAge, …) ─────────────
        matched_scalar = None
        for param in scalar_replacements:
            if re.match(rf"\s*{re.escape(param)}\s*[=~]", stripped):
                matched_scalar = param
                break

        if matched_scalar:
            new_lines.append(
                " " * indent + f"{matched_scalar} = {scalar_replacements[matched_scalar]};\n"
            )
            continue

        new_lines.append(line)

    with open(lphy_path, "w") as fh:
        fh.writelines(new_lines)


def main():
    errors = []

    for r in range(R_START, R_END + 1):
        log_path   = os.path.join(SIM_DATA_DIR, f"sequenceGeneration_r{r}.log")
        trees_path = os.path.join(LPHY_DIR,     f"sequenceGeneration_r{r}_psihealthy.trees")
        lphy_path  = os.path.join(BULK_DIR,     f"bulk_r{r}.lphy")

        if not os.path.exists(lphy_path):
            print(f"[SKIP] {lphy_path} not found")
            continue
        if not os.path.exists(log_path):
            print(f"[SKIP] {log_path} not found")
            continue

        try:
            scalar_replacements: dict[str, str] = {}
            vector_replacements: dict[str, str] = {}

            # ── Scalar params from .log ──────────────────────────────────
            wanted_scalars = {
                p for p, active in TRACE_SCALARS.items()
                if active and p != "rootAge"
            }
            if wanted_scalars:
                scalars = parse_log(log_path, wanted_scalars)
                scalar_replacements.update(scalars)

            # ── rootAge from .trees ──────────────────────────────────────
            if TRACE_SCALARS.get("rootAge"):
                if not os.path.exists(trees_path):
                    print(f"[WARN] r={r}: trees file not found: {trees_path}")
                else:
                    age = get_root_age(trees_path)
                    scalar_replacements["rootAge"] = f"{age:.10g}"

            # ── Vector params from .log ──────────────────────────────────
            for key, cfg in TRACE_VECTORS.items():
                if not cfg["active"]:
                    continue
                vals = parse_vector_from_log(log_path, cfg["prefix"], cfg["size"])
                vector_replacements[cfg["lphy_var"]] = format_vector_lphy(vals)

            # ── Apply to .lphy ───────────────────────────────────────────
            update_lphy(lphy_path, scalar_replacements, vector_replacements)
            summary = {**scalar_replacements,
                       **{k: v[:40] + "...]" for k, v in vector_replacements.items()}}
            print(f"[OK]   r={r}  →  {summary}")

        except Exception as exc:
            msg = f"[ERROR] r={r}: {exc}"
            print(msg, file=sys.stderr)
            errors.append(msg)

    if errors:
        print(f"\n{len(errors)} error(s) encountered.")
        sys.exit(1)
    else:
        print("\nAll replicates processed successfully.")


if __name__ == "__main__":
    main()