#!/bin/bash

#SBATCH -J SM_runMerging
#SBATCH -A uoa04397            # Account code
#SBATCH --time=16:00:00         # Walltime
#SBATCH --cpus-per-task=1       # OpenMP Threads
#SBATCH --ntasks=1              # Not using MPI
#SBATCH --hint=multithread      # A multithreaded job
#SBATCH --mail-user=yil218@aucklanduni.ac.nz
#SBATCH --mail-type=ALL 
chmod +x *sh

./filter.sh
./modifyVCFs.sh
./getGZfiles.sh

module load BCFtools/1.19-GCC-11.3.0
module load Python/3.11.3-gimkl-2022a
python changeVCFnames.py

./mergeVCF.sh
python changePL.py
python getRoorRef.py


