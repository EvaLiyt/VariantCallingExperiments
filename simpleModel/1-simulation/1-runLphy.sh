# Define the path to the slphy executable and the input .lphy file
slphy="/home/yil218/lphy-studio-1.6.1/bin/slphy"
lphy_file="./sequenceGeneration.lphy"
log="./readcount"
result="./simulationData"
nexus="./nexus"
vcf="./vcf"
root="./root"
trees="./trees"
# Create output directories if they don't exist
mkdir -p "${result}" "${nexus}" "${vcf}" "${root}" "${trees}" "${log}"

# Loop to perform the operation 100 times
for i in $(seq 0 99); do
      # Run slphy inside the folder
      ${slphy} -r 1 $lphy_file

      # Rename the generated files to include the folder name (e.g., r1) after the first underscore
      for file in sequenceGeneration_*; do
        if [[ -f "$file" ]]; then
            mv "$file" "$(echo "$file" | sed "s/^sequenceGeneration_/sequenceGeneration_r${i}_/")"
        fi
        done

    # Rename the log file
        for file in sequenceGeneration.log; do
            if [[ -f "$file" ]]; then
                mv "$file" "$(echo "$file" | sed "s/\.log$/_r${i}.log/")"
            fi
        done

     
      mv *.trees ${trees}
      mv *A.nexus ${nexus}
      mv *D.nexus ${root}
      mv *readcount.log ${log}
      mv *.vcf ${vcf}
      
      mv sequenceGeneration_* ${result}
      

      echo "r${i} is done"
done

# Print a message when all simulations are completed
echo "All simulations completed."
