#!/bin/sh
#SBATCH -J CM_runBEAST
#SBATCH --time=8:00:00          # Adjust walltime as needed

# Submit each .sh script
cd xmls
for file in ./*.sh; do
    echo "Submitting $file"
    sbatch "$file"
done
