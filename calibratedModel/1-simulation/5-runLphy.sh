#!/bin/bash
#SBATCH -J CM_Simulation
#SBATCH -A nesi00390            # The account code
#SBATCH --time=120:00:00        # Adjust walltime as needed
#SBATCH --mem=32G              # Adjust memory as needed
#SBATCH --cpus-per-task=4      # Number of CPUs per task
#SBATCH --mail-user=yil218@aucklanduni.ac.nz
#SBATCH --mail-type=ALL

module load Java/17
# Define the path to the slphy executable and the input .lphy file
slphy="/home/yil218/lphy-studio-1.6.1/bin/slphy"
lphy_dir="./lphy_scripts"
result="./simulationData"
log="./readcount"
nexus="./nexus"
vcf="./vcf"
root="./root"
trees="./trees"

# Create output directories if they don't exist
mkdir -p "${log}" "${result}" "${nexus}" "${vcf}" "${root}" "${trees}"

# Loop through each LPhy script
for lphy_file in ${lphy_dir}/*.lphy; do
    # get file name
    basename=$(basename "${lphy_file}")
    id=$(echo "${basename}" | sed -E 's/.*_r([0-9]+)\.lphy/\1/')

    echo "Running simulation for ${id}..."

    # run slphy
    ${slphy} -r 1 "${lphy_file}"
   
    mv ${lphy_dir}/*tree.trees "${trees}/"
    mv ${lphy_dir}/*_A.nexus "${nexus}/"
    mv ${lphy_dir}/*_D.nexus "${root}/"
    mv ${lphy_dir}/*.vcf "${vcf}/"
    mv ${lphy_dir}/*readcount.log "${log}/"
    mv ${lphy_dir}/sequenceGeneration*.trees "${result}/"
    mv ${lphy_dir}/sequenceGeneration*.log "${result}/"

    echo "r${id} simulation completed."
done

echo "All simulations completed."