#!/bin/bash
#SBATCH -J SM_basic_processVCFs
#SBATCH --time=3-00:00:00

chmod +x *sh
module load BCFtools/1.19-GCC-11.3.0
module load Python/3.11.3-gimkl-2022a

./filter.sh
./modifyVCF.sh
./normAndRename.sh

python renameVCF.py ${WORK_DIR:?set WORK_DIR to your own working directory path}/newvcExp/simpleModel/2-variantCalling

./mergeVCFs.sh
python changePL.py
python getRootRef.py

