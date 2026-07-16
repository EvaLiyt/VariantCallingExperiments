#!/bin/bash

#SBATCH -J SM_basic_runMerging
#SBATCH --time=2-00:00:00         # Walltime

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
