#!/bin/bash
#SBATCH -J R_runSimulation
#SBATCH -A nesi00390            # The account code
#SBATCH --time=120:00:00        # Adjust walltime as needed
#SBATCH --mem=32G              # Adjust memory as needed
#SBATCH --cpus-per-task=4      # Number of CPUs per task
#SBATCH --mail-user=yil218@aucklanduni.ac.nz
#SBATCH --mail-type=ALL

module load Java/17
./1-runLphy.sh