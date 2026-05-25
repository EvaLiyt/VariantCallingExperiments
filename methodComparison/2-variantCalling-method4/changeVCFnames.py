import os
import gzip
import shutil
import subprocess

# Define modules to load once, reuse in every command
MODULES = "BCFtools/1.19-GCC-11.3.0 tabix"

def run_with_modules(cmd):
    """Run a shell command with required modules loaded in the same subshell."""
    full_cmd = f"module load {MODULES} && {cmd}"
    subprocess.run(full_cmd, shell=True, check=True, executable='/bin/bash')

# Define input and output directories
input_folder = "./modifiedVCF"
output_folder = "./changedVCF"

# Remove and recreate output folder
if os.path.exists(output_folder):
    shutil.rmtree(output_folder)
    print(f"The directory {output_folder} and all its contents were deleted.")
os.makedirs(output_folder)
print(f"The directory {output_folder} was created.")

# Traverse subfolders in the input folder
for subfolder in os.listdir(input_folder):
    subfolder_path = os.path.join(input_folder, subfolder)
    if os.path.isdir(subfolder_path):
        output_subfolder = os.path.join(output_folder, subfolder)
        os.makedirs(output_subfolder, exist_ok=True)

        for filename in os.listdir(subfolder_path):
            if filename.endswith(".norm.vcf.gz"):
                input_path = os.path.join(subfolder_path, filename)
                output_path = os.path.join(output_subfolder, filename)

                # e.g. "r0_1.sccaller.norm.vcf.gz" -> "r0_1"
                base_name = filename.split(".sccaller")[0]

                compressed_path = output_path + ".gz"

                # Rewrite the VCF header sample names
                with gzip.open(input_path, 'rt') as infile, open(output_path, 'w') as outfile:
                    for line in infile:
                        if line.startswith("#CHROM"):
                            headers = line.strip().split('\t')
                            headers[9:] = [base_name] * len(headers[9:])
                            line = '\t'.join(headers) + '\n'
                        outfile.write(line)

                # Compress, sort, index — each with modules loaded in the same shell
                run_with_modules(f"bgzip -c {output_path} > {compressed_path}")
                os.remove(output_path)
                run_with_modules(f"bcftools sort -Oz -o {compressed_path} {compressed_path}")
                run_with_modules(f"tabix -p vcf {compressed_path}")
                print(f"Processed file: {filename} -> sample name: {base_name}")

print("Finished modifying, compressing, and indexing VCF files.")