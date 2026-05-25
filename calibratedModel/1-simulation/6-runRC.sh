#!/bin/bash
#SBATCH -J CM_runrc
#SBATCH -A uoa04397            # The account code
#SBATCH --time=72:00:00        # Adjust walltime as needed
#SBATCH --mem=2G              # Adjust memory as needed

module load R
Rscript ./6-getReadCounts.R
