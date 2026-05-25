#!/bin/bash
#SBATCH -J CM_processVCFs
#SBATCH -A uoa04397            # Account code
#SBATCH --time=3-00:00:00
#SBATCH --mail-user=yil218@aucklanduni.ac.nz
#SBATCH --mail-type=ALL

chmod +x *sh
module load BCFtools/1.19-GCC-11.3.0
module load Python/3.11.3-gimkl-2022a

./filter.sh
./modifyVCF.sh
./normAndRename.sh
python renameVCF.py /nesi/nobackup/uoa04397/eva/newvcExp/calibratedModel/2-variantCalling
./mergeVCFs.sh
python changePL.py
python getRootRef.py

./6-mapAlignment.sh