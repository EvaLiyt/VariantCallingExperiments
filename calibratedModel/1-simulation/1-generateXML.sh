#!/bin/bash

#SBATCH -J CM_getTrees
#SBATCH -A uoa04397            # Account code
#SBATCH --time=4:00:00         # Walltime
#SBATCH --cpus-per-task=1       # OpenMP Threads
#SBATCH --mem=8G              # Adjust memory as needed
#SBATCH --ntasks=1              # Not using MPI
#SBATCH --hint=multithread      # A multithreaded job
#SBATCH --mail-user=yil218@aucklanduni.ac.nz
#SBATCH --mail-type=ALL

# Load necessary modules
module load Java/17
module load Python/3.11.3-gimkl-2022a
slphy="/home/yil218/lphy-studio-1.6.1/bin/slphy"
lphy_file="./calibrationSimulation.lphy"
##${slphy} $lphy_file

phyloncoVersionFile="/home/yil218/.beast/2.7/phylonco/version.xml"
lphybeastVersionFile="/home/yil218/.beast/2.7/lphybeast/version.xml"
flcVersionFile="/home/yil218/.beast/2.7/lphybeast.flc/version.xml"
lphybeastPath="/home/yil218/beast/bin"

${lphybeastPath}/lphybeast -vf "${phyloncoVersionFile},${lphybeastVersionFile},${flcVersionFile}" -sp -l 20000000 "$lphy_file"

echo "XML files are generated."


