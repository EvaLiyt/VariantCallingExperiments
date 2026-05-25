#!/bin/sh
#SBATCH -J SM_moveFiles
#SBATCH -A uoa04397            # The account code
#SBATCH --time=2:00:00          # Adjust walltime as needed
#SBATCH --mail-user=yil218@aucklanduni.ac.nz
#SBATCH --mail-type=ALL

mkdir -p results
mv xmls/*log results
mv xmls/*trees results
