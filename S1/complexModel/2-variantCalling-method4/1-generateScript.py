import os

base_pileup_dir = "../1-simulation/pileup"
vcf_base = "../1-simulation/vcf"
bulk_dir = "../2-variantCalling/bulkSam"
jobs_dir = "jobs"
fasta = os.path.expanduser("~/data/10K.fna")
sccaller = "../../sccaller_labVersion.py"
out_base = "step1"

os.makedirs(jobs_dir, exist_ok=True)


for repeat in range(100):
    repeat_name = f"r{repeat}"
    repeat_dir = os.path.join(base_pileup_dir, repeat_name)
    out_dir = os.path.join(out_base, repeat_name)
    cell_names = [
        os.path.splitext(f)[0]
        for f in sorted(os.listdir(repeat_dir))
        if f.endswith(".pileup")
    ]
    cmds = []
    for cell in cell_names:
        pileup_file = os.path.join(repeat_dir, f"{cell}.pileup")
        snp_file = os.path.join(vcf_base, f"sequenceGeneration_{repeat_name}_snp.vcf")
        bulk_file = os.path.join(bulk_dir, f"{repeat_name}.sorted.bam")
        output_file = os.path.join(out_dir, f"{repeat_name}_{cell}.sccaller.vcf")
        cmd = f"python {sccaller} --pileup {pileup_file} --fasta {fasta} --snp_in {snp_file} --output {output_file} --snp_type hsnp --bulk {bulk_file} --min_depth 2 --minvar 0 --bias 0.6 --lamb 2000 --mapq 30 --engine samtools --cpu_num 8"
        cmds.append(cmd)

    lines = []
    lines.append("#!/bin/bash")
    lines.append(f"#SBATCH -J step1_RM_{repeat_name}")
    lines.append("#SBATCH --time=3:00:00")
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
    print(f"Script written: {script_path}")
    os.system(f"sbatch {script_path}") 

