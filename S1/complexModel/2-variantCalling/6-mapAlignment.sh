#!/bin/bash

#SBATCH -J RM_mapAlignment
#SBATCH --time=3-00:00:00

chmod +x ./getMappedFile.sh
module load Python

./getMappedFile.sh
python replaceQuestionMarks.py

