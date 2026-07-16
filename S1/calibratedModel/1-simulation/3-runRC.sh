#!/bin/bash
#SBATCH -J SM_runrc
#SBATCH --time=72:00:00        # Adjust walltime as needed
#SBATCH --mem=4G              # Adjust memory as needed
#SBATCH --cpus-per-task=1
module load R
Rscript ./3-getReadCounts.R
