#!/bin/bash
#SBATCH -J SM_replace?s
#SBATCH -A nesi00390            # The account code
#SBATCH --time=2:00:00        # Adjust walltime as needed
#SBATCH --mem=16G              # Adjust memory as needed
#SBATCH --cpus-per-task=4      # Number of CPUs per task
#SBATCH --mail-user=yil218@aucklanduni.ac.nz
#SBATCH --mail-type=ALL

module load Python/3.11.3-gimkl-2022a
python 6-replaceQuestionMarks.py
