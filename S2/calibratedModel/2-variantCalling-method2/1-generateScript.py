#!/usr/bin/env python3
"""
Generate SLURM submission scripts for sccaller variant calling across multiple repeats.

Expects directory structure:
  ../1-simulation/pileup/r<N>/*.pileup
  ../1-simulation/vcf/sequenceGeneration_r<N>_snp.vcf

Outputs all scripts to current directory:
  ./run_r<N>.sh
"""

import os
import sys
import textwrap
from pathlib import Path

# Configuration
NUM_REPEATS = 100
BASE_PILEUP_DIR = "../1-simulation/pileup"
VCF_BASE_DIR = "../1-simulation/vcf"
REFERENCE = "~/data/10K.fna"

# SLURM settings
SLURM_TIME = "2:00:00"
SLURM_MEM = "4G"

# Sccaller settings
SCCALLER_SCRIPT = "../../sccaller_labVersion.py"
MIN_DEPTH = 1
MIN_VAR = 1
ENGINE = "samtools"
TARGET_TYPE = "hsnp"

def validate_paths():
    """Validate that base directories exist."""
    if not os.path.isdir(BASE_PILEUP_DIR):
        print(f"[ERROR] Base pileup directory not found: {BASE_PILEUP_DIR}", file=sys.stderr)
        sys.exit(1)
    if not os.path.isdir(VCF_BASE_DIR):
        print(f"[ERROR] VCF base directory not found: {VCF_BASE_DIR}", file=sys.stderr)
        sys.exit(1)
    print(f"[OK] Directories validated")

def get_repeat_dirs():
    """
    Find all r<N> directories in the pileup base.
    Returns sorted list of (repeat_num, repeat_dir, repeat_name) tuples.
    """
    repeats = []
    try:
        for item in os.listdir(BASE_PILEUP_DIR):
            if not item.startswith("r"):
                continue
            try:
                repeat_num = int(item[1:])
                repeat_dir = os.path.join(BASE_PILEUP_DIR, item)
                if os.path.isdir(repeat_dir):
                    repeats.append((repeat_num, repeat_dir, item))
            except ValueError:
                # Not a valid r<N> format
                continue
    except OSError as e:
        print(f"[ERROR] Failed to list {BASE_PILEUP_DIR}: {e}", file=sys.stderr)
        sys.exit(1)

    repeats.sort(key=lambda x: x[0])
    return repeats

def find_pileup_files(repeat_dir):
    """Find all .pileup files in the repeat directory."""
    try:
        files = [f for f in os.listdir(repeat_dir) if f.endswith(".pileup")]
        return sorted(files)
    except OSError as e:
        print(f"[WARN] Failed to list {repeat_dir}: {e}", file=sys.stderr)
        return []

def verify_input_files(repeat_dir, repeat_name, pileup_files):
    """Verify that required input files exist."""
    vcf_path = os.path.join(VCF_BASE_DIR, f"sequenceGeneration_{repeat_name}_snp.vcf")

    issues = []

    if not pileup_files:
        issues.append(f"  - No .pileup files found in {repeat_dir}")

    if not os.path.exists(vcf_path):
        issues.append(f"  - VCF file not found: {vcf_path}")

    if not os.path.exists(os.path.expanduser(REFERENCE)):
        issues.append(f"  - Reference file not found: {REFERENCE}")

    if not os.path.exists(SCCALLER_SCRIPT):
        issues.append(f"  - sccaller script not found: {SCCALLER_SCRIPT}")

    return issues

def generate_slurm_header(repeat_name):
    """Generate SLURM header for the submission script."""
    return textwrap.dedent(f"""\
        #!/bin/bash
        #SBATCH -J CM_{repeat_name}
        #SBATCH --time={SLURM_TIME}
        #SBATCH --mem={SLURM_MEM}
        #SBATCH -o {repeat_name}_%j.out
        #SBATCH -e {repeat_name}_%j.err
        """)

def generate_module_loads():
    """Generate module load commands."""
    return [
        "module load Python/3.11.3-gimkl-2022a",
        "module load SAMtools",
        ""
    ]

def generate_sccaller_command(pileup_path, repeat_name, cell_name, vcf_path, output_dir="vsResults"):
    """Generate sccaller invocation for a single pileup file."""
    # Output to csResults/{repeat_name}/{repeat_name}_{cell_name}.vcf
    command = (
        f"python3 {SCCALLER_SCRIPT} "
        f"--min_depth {MIN_DEPTH} "
        f"--minvar {MIN_VAR} "
        f"--engine {ENGINE} "
        f"-p {pileup_path} "
        f"-f {REFERENCE} "
        f"-s {vcf_path} "
        f"-t {TARGET_TYPE} "
        f"-o {output_dir}/{repeat_name}/{repeat_name}_{cell_name}.vcf"
    )
    return command

def generate_script(repeat_name, repeat_dir, pileup_files, vcf_path, output_dir="csResults"):
    """
    Generate a complete SLURM submission script for a single repeat.
    Returns the script content as a string.
    """
    lines = []

    # Header
    lines.append(generate_slurm_header(repeat_name))

    # Comments
    lines.append(f"# Variant calling for {repeat_name}")
    lines.append(f"# Processing {len(pileup_files)} pileup file(s)")
    lines.append("")

    # Module loads
    lines.extend(generate_module_loads())

    # Set up error handling
    lines.append("set -e  # Exit on first error")
    lines.append("set -u  # Exit if undefined variable")
    lines.append("")

    # Create output directories
    lines.append(f"mkdir -p {output_dir}/{repeat_name}")
    lines.append("")

    # Generate commands for each pileup file
    lines.append("# Variant calling commands")
    for pileup_file in pileup_files:
        cell_name = os.path.splitext(pileup_file)[0]
        pileup_path = os.path.join(repeat_dir, pileup_file)
        command = generate_sccaller_command(pileup_path, repeat_name, cell_name, vcf_path, output_dir)
        lines.append("")
        lines.append(f"echo '[{repeat_name}] Processing {cell_name}...'")
        lines.append(command)

    lines.append("")
    lines.append(f"echo '[{repeat_name}] Complete'")

    return "\n".join(lines)

def main():
    """Main entry point."""
    print("Generating SLURM submission scripts for sccaller variant calling\n")

    # Validate base directories
    validate_paths()

    # Find all repeat directories
    repeats = get_repeat_dirs()

    if not repeats:
        print("[WARN] No r<N> directories found", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(repeats)} repeat(s): {', '.join(r[2] for r in repeats)}\n")

    # Create scripts output folder
    scripts_dir = "scripts"
    output_dir = "csResults"
    os.makedirs(scripts_dir, exist_ok=True)
    print(f"Output scripts folder: ./{scripts_dir}/")
    print(f"Output results folder: ./{output_dir}/<repeat_name>/\n")

    success_count = 0
    skip_count = 0

    for repeat_num, repeat_dir, repeat_name in repeats:
        print(f"Processing {repeat_name}...")

        # Find pileup files
        pileup_files = find_pileup_files(repeat_dir)

        # Build VCF path
        vcf_path = os.path.join(VCF_BASE_DIR, f"sequenceGeneration_{repeat_name}_snp.vcf")

        # Verify inputs
        issues = verify_input_files(repeat_dir, repeat_name, pileup_files)
        if issues:
            print(f"  [WARN] Skipping {repeat_name} due to missing inputs:")
            for issue in issues:
                print(issue, file=sys.stderr)
            skip_count += 1
            continue

        # Generate script
        script_content = generate_script(repeat_name, repeat_dir, pileup_files, vcf_path, output_dir)

        # Write script to scripts folder
        script_path = os.path.join(scripts_dir, f"run_{repeat_name}.sh")
        try:
            with open(script_path, "w") as f:
                f.write(script_content)
            # Make executable
            os.chmod(script_path, 0o755)
            print(f"  [OK] {script_path} ({len(pileup_files)} pileup file(s))")
            success_count += 1
        except OSError as e:
            print(f"  [ERROR] Failed to write {script_path}: {e}", file=sys.stderr)
            skip_count += 1

    print(f"\n{'='*60}")
    print(f"Summary: {success_count} script(s) generated, {skip_count} skipped")
    print(f"Scripts location: ./{scripts_dir}/")
    print(f"{'='*60}")

    if success_count == 0:
        print("[ERROR] No scripts were generated", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
