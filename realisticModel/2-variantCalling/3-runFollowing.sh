#!/bin/bash

#SBATCH -J RM_filter&Modify
#SBATCH -A nesi00390            # Account code
#SBATCH --time=16:00:00         # Walltime
#SBATCH --cpus-per-task=1       # OpenMP Threads
#SBATCH --ntasks=1              # Not using MPI
#SBATCH --hint=multithread      # A multithreaded job
#SBATCH --mail-user=yil218@aucklanduni.ac.nz
#SBATCH --mail-type=ALL 

./filter.sh
./modifyVCFs.sh



