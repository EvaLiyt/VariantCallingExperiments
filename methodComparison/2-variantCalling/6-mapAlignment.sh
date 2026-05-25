#!/bin/bash

#SBATCH -J CM_mapAlignment
#SBATCH -A uoa04397            # Account code
#SBATCH --time=3-00:00:00

chmod +x ./getMappedFile.sh
module load Python

./getMappedFile.sh
python replaceQuestionMarks.py

