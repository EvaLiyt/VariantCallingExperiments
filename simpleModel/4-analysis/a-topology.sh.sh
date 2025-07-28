#!/bin/bash
#SBATCH -J runTopologyCheck
#SBATCH -A uoa04397
#SBATCH --time=5:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=yil218@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH -o slurm_r${id}.out
#SBATCH -e slurm_r${id}.err

module load Java/17
chmod +x a-checkTopology.sh
./a-checkTopology.sh
