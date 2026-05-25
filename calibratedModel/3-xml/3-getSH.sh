#!/bin/sh
#SBATCH -J CM_getSHs
#SBATCH -A uoa04397            # The account code

# Define the list of folders
folders=("ErrModelScripts" "wtErrModelScripts" "TrueErrModelScripts" "TruewtErrModelScripts")

# Loop over each folder
for folder in "${folders[@]}"
do
    # Extract the xml_type from folder name (remove "Scripts")
    xml_type=$(basename "$folder" | sed 's/Scripts//')

    # Loop over r0 to r99
    for r in {0..99}
    do
        # Output .sh file path
        sh_file="xmls/${xml_type}_r${r}.sh"

        # Create the .sh script
        echo "#!/bin/sh" > "$sh_file"
        echo "" >> "$sh_file"
        echo "#SBATCH -J CM_${xml_type}_r${r}      # The job name" >> "$sh_file"
        echo "#SBATCH -A uoa04397                # The account code" >> "$sh_file"
        echo "#SBATCH --time=5-00:00:00             # The walltime" >> "$sh_file"
        echo "#SBATCH --mem=3G                    # Total memory" >> "$sh_file"
        echo "#SBATCH --cpus-per-task=2           # OpenMP Threads" >> "$sh_file"
        echo "#SBATCH --ntasks=1                  # Not using MPI" >> "$sh_file"
        echo "#SBATCH -D ./                       # Initial directory" >> "$sh_file"
        echo "#SBATCH -o out/${xml_type}_r${r}.out    # Output file" >> "$sh_file"
        echo "#SBATCH -e out/${xml_type}_r${r}.err    # Error file" >> "$sh_file"

        echo "" >> "$sh_file"
        echo "# sacct -j JOBID --format=\"ReqMem,MaxRSS,CPUTime,AveCPU,Elapsed\"" >> "$sh_file"
        echo "module load beagle-lib/4.0.0-GCC-11.3.0" >> "$sh_file"
        echo "module load Java/17" >> "$sh_file"
        echo "" >> "$sh_file"

        # Define XML filename
        xml_file="./${xml_type}_r${r}.xml"
        echo "# Running beast on $xml_file" >> "$sh_file"
        echo "srun /home/yil218/beast-2.7/bin/beast -beagle_SSE $xml_file" >> "$sh_file"
    
    done
done
