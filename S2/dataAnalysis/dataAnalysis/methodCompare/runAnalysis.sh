#!/bin/bash
#SBATCH -J runAnalysis
#SBATCH --time=24:00:00
#SBATCH --mem=2G

module load Python
python analyseData.py
python plot.py