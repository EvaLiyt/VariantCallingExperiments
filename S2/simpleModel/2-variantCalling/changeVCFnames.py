import os
import gzip
import shutil
import subprocess

try:
    subprocess.run("module load BCFtools/1.19-GCC-11.3.0", shell=True, check=True, executable='/bin/bash')
    print("Successfully loaded bcftools module")
except subprocess.CalledProcessError as e:
    print(f"Error loading bcftools module: {e}")
    
try:
    subprocess.run("module load tabix", shell=True, check=True, executable='/bin/bash')
    print("Successfully loaded tabix module")
except subprocess.CalledProcessError as e:
    print(f"Error loading tabix module: {e}")

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
        # Create corresponding subfolder in the output folder
        output_subfolder = os.path.join(output_folder, subfolder)
        os.makedirs(output_subfolder, exist_ok=True)

        # Traverse all .vcf.gz.norm.vcf.gz files in the subfolder
        for filename in os.listdir(subfolder_path):
            if filename.endswith(".norm.vcf.gz"):
                input_path = os.path.join(subfolder_path, filename)
                output_path = os.path.join(output_subfolder, filename)

                # Extract base filename without extension
                base_name = filename.replace(".norm.vcf.gz", "")

                # Define paths for bgzip and tabix
                compressed_path = output_path + ".gz"
                index_path = compressed_path + ".tbi"

                # Process the VCF file
                with gzip.open(input_path, 'rt') as infile, open(output_path, 'w') as outfile:
                    for line in infile:
                        if line.startswith("#CHROM"):
                            headers = line.strip().split('\t')
                            # Replace sample columns starting from the 10th with base_name
                            headers[9:] = [f"{base_name}" for _ in range(len(headers[9:]))]
                            line = '\t'.join(headers) + '\n'
                        outfile.write(line)

                # Compress the output file using bgzip
                bgzip_command = f"bgzip -c {output_path} > {compressed_path}"
                subprocess.run(bgzip_command, shell=True, check=True)

                # Remove the uncompressed file
                os.remove(output_path)
                sort_command = f"bcftools sort -Oz -o {compressed_path} {compressed_path}"
                subprocess.run(sort_command, shell=True, check=True)
				
                # Generate the index file using tabix
                tabix_command = f"tabix -p vcf {compressed_path}"
                subprocess.run(tabix_command, shell=True, check=True)

                print(f"Processed file: {filename}")

print("Finished modifying, compressing, and indexing VCF files.")
