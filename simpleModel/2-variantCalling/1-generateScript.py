import os
import random
import textwrap

# Set seed for reproducibility if needed
random.seed(42)

# Number of repeats
num_repeats = 100

# Base path where folders are located
base_pileup_dir = "../1-simulation/pileup"
vcf_base = "../1-simulation/vcf"
reference = "~/data/chr1.fna"


for r in range(100):
    repeat_name = f"r{r}"
    repeat_dir = os.path.join(base_pileup_dir, repeat_name)

    if not os.path.isdir(repeat_dir):
        print(f"Warning: Directory {repeat_dir} not found. Skipping.")
        continue

    # Randomly choose parameters
    min_depth = random.randint(2, 8)
    min_var = random.randint(1, min_depth - 1)

   # Write parameter file
    parameter_file_path = os.path.join(repeat_dir, f"parameters_{repeat_name}.txt")
    with open(parameter_file_path, "w") as param_file:
        param_file.write(f"min_depth: {min_depth}\n")
        param_file.write(f"min_var: {min_var}\n")

    slurm_header = textwrap.dedent(f"""\
    #!/bin/bash
    #SBATCH -J SM_{repeat_name}
    #SBATCH -A nesi00390
    #SBATCH --time=72:00:00
    #SBATCH --mem=8G
    #SBATCH --cpus-per-task=4
    """)
    # Add email options only for r99
    if repeat_name == "r99":
        slurm_header += textwrap.dedent("""\
#SBATCH --mail-user=yil218@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
""")


    bash_lines = [slurm_header, f"# Bash script for {repeat_name}", f"module load Python/3.11.3-gimkl-2022a",  f"# min_depth = {min_depth}, min_var = {min_var}\n"]

    # Find all .pileup files in the directory
    pileup_files = [f for f in os.listdir(repeat_dir) if f.endswith(".pileup")]

    for pileup_file in pileup_files:
        cell_name = os.path.splitext(pileup_file)[0]  # remove '.pileup'
        pileup_path = os.path.join(repeat_dir, pileup_file)
        output_vcf = f"{repeat_dir}/{repeat_name}_{cell_name}.vcf"
        vcf_path = os.path.join(vcf_base, f"sequenceGeneration_{repeat_name}_snp.vcf")

        command = (
            f"python3 ../../sccaller_pileup.py "
            f"--min_depth {min_depth} --minvar {min_var} "
            f"-p {pileup_path} "
            f"-f {reference} "
            f"-s {vcf_path} "
            f"-t hsnp "
            f"-o {output_vcf}"
        )
        bash_lines.append(command)

    # Write the bash script in the pileup directory
    script_path = os.path.join(repeat_dir, f"run_{repeat_name}.sh")
    with open(script_path, "w") as f:
        f.write("\n".join(bash_lines))

    print(f"Script written: {script_path}")