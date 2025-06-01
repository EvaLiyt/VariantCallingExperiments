#!/bin/bash
#SBATCH -J RM_runGZ
#SBATCH -A nesi00390            # Account code
#SBATCH --time=16:00:00         # Walltime
#SBATCH --cpus-per-task=1       # OpenMP Threads
#SBATCH --ntasks=1              # Not using MPI
#SBATCH --hint=multithread      # A multithreaded job
#SBATCH --mail-user=yil218@aucklanduni.ac.nz
#SBATCH --mail-type=ALL 

# Define input and reference
input_folder="./modifiedVCF"
reference_fasta="~/data/chr1.fna"

# Loop through each rX folder
for subfolder in "${input_folder}"/*; do
    if [ -d "$subfolder" ]; then
        subfolder_name=$(basename "$subfolder")
        slurm_script="${subfolder}/bcftools_${subfolder_name}.slurm"

        cat > "$slurm_script" <<EOL
#!/bin/bash
#SBATCH -J getgz_${subfolder_name}
#SBATCH -A nesi00390
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --mail-user=${mail}
#SBATCH --mail-type=END,FAIL
#SBATCH -o ${subfolder}/bcftools_${subfolder_name}.out
#SBATCH -e ${subfolder}/bcftools_${subfolder_name}.err

module load GCCcore/11.3.0
module load BCFtools/1.19-GCC-11.3.0

echo "Processing ${subfolder_name}"

for file in "${subfolder}"/*.vcf; do
    filename=\$(basename "\${file}")
    gz_file="\${subfolder}/\${filename%.vcf}.vcf.gz"
    norm_file="\${subfolder}/\${filename%.vcf}.norm.vcf.gz"

    bcftools view "\${file}" -Oz -o "\${gz_file}"
    bcftools norm -Oz -f ${reference_fasta} -o "\${norm_file}" "\${gz_file}"
    bcftools index "\${norm_file}"
done

echo "Finished ${subfolder_name}"
EOL

        echo "Generated: $slurm_script"

        sbatch "$slurm_script"
    fi
done

echo "All SLURM scripts generated."
