#!/bin/bash

#SBATCH -J SM_mapAlignment
#SBATCH --time=3-00:00:00
#SBATCH --hint=multithread      # A multithreaded job

chmod +x ./getMappedFile.sh
module load Python

./getMappedFile.sh
python replaceQuestionMarks.py

