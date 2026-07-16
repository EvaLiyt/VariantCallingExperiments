#!/bin/sh
#SBATCH -J SM_moveFiles
#SBATCH --time=2:00:00          # Adjust walltime as needed

mkdir -p results
mv xmls/*log results
mv xmls/*trees results
