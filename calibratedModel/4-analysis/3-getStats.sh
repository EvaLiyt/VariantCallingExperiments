#!/bin/bash

#SBATCH -J CM_getStats
#SBATCH -A uoa04397            # The account code
#SBATCH --time=12:00:00        # Adjust walltime as needed
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2	# OpenMP Threads
#SBATCH --ntasks=1              # not use MPI
#SBATCH --hint=multithread	# A multithreaded job

module load Python/3.11.3-gimkl-2022a
# Run the Python script
python ./2-stats.py
