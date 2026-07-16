#!/bin/bash
# Array of minD values
minD=(1 2 3 4 5)
# Create folders DxVx where D = each value in minD, V = 1..D
for D in "${minD[@]}"; do
    for ((V=1; V<=D; V++)); do
        folder="D${D}V${V}"
        mkdir -p "$folder/jobs"
        cat > "$folder/1-generateScript.py" << PYEOF
import os
# --- Parameters for this folder ---
min_depth = ${D}
min_var   = ${V}
# --- Paths (relative to this DxVx folder) ---
base_pileup_dir = "../../1-simulation/pileup"
vcf_base        = "../../1-simulation/vcf"
bulk_dir        = "../bulkSam"
jobs_dir        = "jobs"
fasta           = "$HOME/data/10K.fna"
sccaller        = "../../../sccaller_labVersion.py"
out_base        = "vcResults"
os.makedirs(jobs_dir, exist_ok=True)
for repeat in range(100):
    repeat_name = f"r{repeat}"
    repeat_dir  = os.path.join(base_pileup_dir, repeat_name)
    out_dir     = os.path.join(out_base, repeat_name)
    # Discover cell names from pileup filenames
    cell_names = [
        os.path.splitext(f)[0]
        for f in sorted(os.listdir(repeat_dir))
        if f.endswith(".pileup")
    ]
    cmds = []
    for cell in cell_names:
        pileup_file = os.path.join(repeat_dir, f"{cell}.pileup")
        snp_file    = os.path.join(vcf_base, f"sequenceGeneration_{repeat_name}_snp.vcf")
        bulk_file   = os.path.join(bulk_dir, f"{repeat_name}.sorted.bam")
        output_file = os.path.join(out_dir, f"{repeat_name}_{cell}.vcf")
        log_file    = os.path.join(out_dir, f"sc_{repeat_name}_{cell}_D{min_depth}V{min_var}to-1.log")
        cmd = f"python {sccaller} --pileup {pileup_file} --fasta {fasta} --snp_in {snp_file} --output {output_file} --snp_type hsnp --bulk {bulk_file} --bulk_min_var 5 --min_depth {min_depth} --minvar {min_var} --bias 0.75 --engine samtools --cpu_num 2 --format vcf > {log_file} 2>&1 &"
        cmds.append(cmd)
    lines = []
    lines.append("#!/bin/bash")
    lines.append(f"#SBATCH -J D{min_depth}V{min_var}_RM_{repeat_name}")
    lines.append("#SBATCH --time=6:00:00")
    lines.append("#SBATCH --mem=4G")
    lines.append("")
    lines.append("module load Python")
    lines.append("module load SAMtools")
    lines.append(f"mkdir -p {out_dir}")
    lines.append("")
    lines.extend(cmds)
    lines.append("")
    lines.append("wait")
    lines.append("")
    script_content = "\n".join(lines)
    script_path = os.path.join(jobs_dir, f"run_{repeat_name}.sh")
    with open(script_path, "w") as f:
        f.write(script_content)
    print(f"Script written: {script_path} ({len(cell_names)} cells)")
PYEOF
        echo "Created: $folder/1-generateScript.py"
    done
done
echo ""
echo "Done. Folders created:"
for D in "${minD[@]}"; do
    for ((V=1; V<=D; V++)); do
        echo "  D${D}V${V}/"
    done
done