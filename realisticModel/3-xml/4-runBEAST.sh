#!/bin/sh
#SBATCH -J RM_runBEAST
#SBATCH -A nesi00390            # The account code
#SBATCH --time=8:00:00          # Adjust walltime as needed
#SBATCH --mem=8G                # Adjust memory as needed
#SBATCH --cpus-per-task=4       # Number of CPUs per task
#SBATCH --mail-user=yil218@aucklanduni.ac.nz
#SBATCH --mail-type=ALL

chmod +x ./xmls/*.sh

# Submit each .sh script
cd xmls
for file in ./*.sh; do
    echo "Submitting $file"
    sbatch "$file"
done
