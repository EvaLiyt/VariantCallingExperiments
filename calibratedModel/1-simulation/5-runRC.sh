#!/bin/bash
#SBATCH -J CM_runrc
#SBATCH -A uoa04397            # The account code
#SBATCH --time=24:00:00        # Adjust walltime as needed
#SBATCH --mem=8G              # Adjust memory as needed
#SBATCH --cpus-per-task=4      # Number of CPUs per task
#SBATCH --mail-user=yil218@aucklanduni.ac.nz
#SBATCH --mail-type=ALL

module load R
Rscript ./5-getReadCounts.R
