#!/bin/bash

#SBATCH -J SM_mapAlignment
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=1       # OpenMP Threads
#SBATCH --ntasks=1              # Not using MPI
#SBATCH --hint=multithread      # A multithreaded job

chmod +x ./getMappedFile.sh
module load Python

./getMappedFile.sh
python replaceQuestionMarks.py

