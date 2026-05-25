#!/bin/bash
#SBATCH -J runAnalysis
#SBATCH -A uoa04397
#SBATCH --time=24:00:00
#SBATCH --mem=2G
#SBATCH --mail-user=yil218@aucklanduni.ac.nz
#SBATCH --mail-type=ALL

module load Python
python analyseData.py
python plot.py