#!/bin/bash
#SBATCH -J SM_runrc
#SBATCH -A uoa04397            # The account code
#SBATCH --time=1:00:00        # Adjust walltime as needed

module load R
Rscript ./2-getReadCounts.R
