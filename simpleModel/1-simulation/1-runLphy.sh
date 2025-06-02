#!/bin/bash

# Paths
lphy_dir="./lphy_scripts"
output_dir="./slurm_jobs"
mkdir -p "${output_dir}"

# Loop through each LPhy script and create a corresponding SLURM job script
for lphy_file in ${lphy_dir}/*.lphy; do
    basename=$(basename "${lphy_file}")
    id=$(echo "${basename}" | sed -E 's/.*_r([0-9]+)\.lphy/\1/')

    job_script="${output_dir}/run_r${id}.sh"

    cat > "${job_script}" <<EOF
#!/bin/bash
#SBATCH -J SM_r${id}
#SBATCH -A nesi00390
#SBATCH --time=5:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=yil218@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH -o slurm_r${id}.out
#SBATCH -e slurm_r${id}.err

module load Java/17

slphy="/home/yil218/lphy-studio-1.6.1/bin/slphy"
lphy_file="${lphy_file}"

result="./simulationData"
log="./readcount"
nexus="./nexus"
vcf="./vcf"
root="./root"
trees="./trees"

mkdir -p "\${log}" "\${result}" "\${nexus}" "\${vcf}" "\${root}" "\${trees}"

echo "Running simulation for r${id}..."
\${slphy} -r 1 "\${lphy_file}"

mv "\${lphy_dir}"/*.trees "\${trees}/"
mv "\${lphy_dir}"/*A.nexus "\${nexus}/"
mv "\${lphy_dir}"/*D.nexus "\${root}/"
mv "\${lphy_dir}"/*.vcf "\${vcf}/"
mv "\${lphy_dir}"/*readcount.log "\${log}/"
mv "\${lphy_dir}"/sequenceGeneration*.log "\${result}/"

echo "r${id} simulation completed."
EOF

    # Submit the job
    sbatch "${job_script}"
    echo "Submitted job for r${id}"
done
