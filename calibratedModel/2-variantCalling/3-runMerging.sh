#!/bin/bash

#SBATCH -J CM_runMerging
#SBATCH -A uoa04397            # Account code
#SBATCH --time=72:00:00         # Walltime
#SBATCH --mail-user=yil218@aucklanduni.ac.nz
#SBATCH --mail-type=ALL 
chmod +x *sh
rm *out
rm *log
./filter.sh
./modifyVCFs.sh
./getGZfiles.sh

module load BCFtools/1.19-GCC-11.3.0
module load Python/3.11.3-gimkl-2022a
python changeVCFnames.py

./mergeVCF.sh
python changePL.py
python getRootRef.py

./4-generateMappedFile.sh
python 5-replaceQuestionMarks.py
