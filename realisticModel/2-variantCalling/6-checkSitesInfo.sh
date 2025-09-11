#!/bin/bash
#SBATCH -J RM_testInfoSites
#SBATCH -A uoa04397
#SBATCH --time=5:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=4

slphy="/home/yil218/lphy-studio-1.6.1/bin/slphy"
summary_file="summary_infoSites.txt"
for r in {0..99}
do
  cat > script_r${r}.lphy <<EOF
data {
    // alignment
    RA = readFasta(file = "alignment/alignment_r${r}.fasta", sequenceType=phasedGenotype());
    var = variableSites(alignment=RA);
    inf = informativeSites(alignment=RA);
}
EOF

  # Run and append filtered output to summary file
  echo "=== r${r} ===" >> "$summary_file"
  ${slphy} script_r${r}.lphy | grep -E "variable sites|informative sites" >> "$summary_file"
done
