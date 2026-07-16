#!/bin/bash
#SBATCH -J RM_runrc
#SBATCH --time=72:00:00        # Adjust walltime as needed
#SBATCH --mem=2G              # Adjust memory as needed

module load R
Rscript ./3-getReadCounts.R
