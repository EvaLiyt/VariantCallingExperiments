#!/bin/bash
#SBATCH -J runART
#SBATCH -A uoa04397
#SBATCH --time=24:00:00
#SBATCH --mem=4G

art=/home/yil218/art_bin_MountRainier/art_illumina

mkdir -p bulkSeqResult

# Run ART on each FASTA at 50x
for fasta_file in ../1-simulation/root/*.fasta; do

    id=$(basename "$fasta_file" .fasta | sed 's/sequenceGeneration_\(.*\)_FD/\1/')

    $art -ss HS25 -i "$fasta_file" -o "bulkSeqResult/${id}_pairend" -l 150 -f 25 -p -m 500 -s 30 -sam
done
