#!/bin/bash
#SBATCH -J processVCFs
#SBATCH --time=3-00:00:00

chmod +x *sh
module load BCFtools/1.19-GCC-11.3.0
module load Python/3.11.3-gimkl-2022a

./filter.sh
./modifyVCF.sh
./normAndRename.sh

./mergeVCFs.sh
python changePL.py
python getRootRef.py

./6-mapAlignment.sh