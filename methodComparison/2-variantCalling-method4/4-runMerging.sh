#!/bin/bash

#SBATCH -J CM_M4_runMerging
#SBATCH -A uoa04397            # Account code
#SBATCH --time=2-00:00:00         # Walltime
#SBATCH --mail-user=yil218@aucklanduni.ac.nz
#SBATCH --mail-type=ALL 
chmod +x *sh


./filter.sh
./modifyVCFs.sh
./getGZfiles.sh

module load BCFtools/1.19-GCC-11.3.0
module load Python
python changeVCFnames.py

./mergeVCF.sh
python changePL.py
python getQuestionMarks.py

./5-generateMappedFile.sh