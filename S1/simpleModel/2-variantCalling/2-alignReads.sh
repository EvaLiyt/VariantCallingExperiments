#!/bin/bash
#SBATCH -J indexSam
#SBATCH --time=72:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2

chmod +x changeSamName.sh
./changeSamName.sh

mkdir -p ./bulkSam
module load SAMtools/1.16.1-GCC-11.3.0

for x in $(seq 0 99); do
    # Convert SAM to BAM, sort, and index in one go
    samtools sort -@ 2 -o bulkSam/r${x}.sorted.bam bulkSeqResult/r${x}_pairend.sam
    samtools index bulkSam/r${x}.sorted.bam
done
