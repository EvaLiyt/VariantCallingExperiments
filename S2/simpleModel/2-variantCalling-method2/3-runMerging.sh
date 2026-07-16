#!/bin/bash

#SBATCH -J SM_runMerging
#SBATCH --time=72:00:00         # Walltime
chmod +x *sh
rm *out
rm *err
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

python mergeVCF.py

./4-generateMappedFile.sh
python 5-replaceQuestionMarks.py

