#!/bin/bash

#SBATCH -J runParameterGetting
#SBATCH -A nesi00390            # The account code
#SBATCH --time=32:00:00        # Adjust walltime as needed
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2	# OpenMP Threads
#SBATCH --ntasks=1              # not use MPI
#SBATCH --hint=multithread	# A multithreaded job
#SBATCH --mail-user=yil218@aucklanduni.ac.nz
#SBATCH --mail-type=ALL

module load R
# Run the Python script
Rscript ./7-extractParameters.R