#!/bin/bash

#SBATCH -J RM_getConsensusTrees
#SBATCH -A uoa04397            # The account code
#SBATCH --time=12:00:00        # Adjust walltime as needed
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2	# OpenMP Threads
#SBATCH --ntasks=1              # not use MPI
#SBATCH --hint=multithread	# A multithreaded job
#SBATCH --mail-user=yil218@aucklanduni.ac.nz
#SBATCH --mail-type=ALL


# Run the Python script
cd /home/yil218/beast-2.7/bin
folder="/nesi/nobackup/uoa04397/eva/vcExp/realisticModel/3-xml/results"

for filename in ${folder}/*trees; do
    echo "processing ${filename}"
  ./treeannotator -topology CCD0 -height CA ${filename} "${filename}.txt"
  echo "done ${filename}"
done
