#!/bin/bash
#SBATCH -J CM_m2_mapAlignment
#SBATCH -A uoa04397
#SBATCH --time=3:00:00

module load Python

# copy root_ref fasta from each D*V* (excluding D1V1) into ./D*V*/merge/
for folder in ../2-variantCalling/D*V*/; do
    fname=$(basename "$folder")
    [[ "$fname" == "D1V1" ]] && continue

    mkdir -p "./${fname}/merge"
    cp "${folder}/merge/root_ref"*.fasta "./${fname}/merge/"
done

chmod +x ./getMappedFile.sh
./getMappedFile.sh

python replaceQuestionMarks.py
