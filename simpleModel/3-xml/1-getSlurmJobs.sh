#!/bin/bash
#SBATCH -J SM_runLphybeast
#SBATCH -A uoa04397
#SBATCH --time=5:00:00

# Define lphybeast executable and version files
lphybeastPath="/home/yil218/beast-2.7/bin"

# Define your target directories
directories=( "./ErrModelScripts" "./wtErrModelScripts" "./TrueErrModelScripts" "./TruewtErrModelScripts" )

# Loop through each directory
for dir in "${directories[@]}"; do
    for lphy_file in "${dir}"/*.lphy; do
        [ -e "$lphy_file" ] || continue  # skip if no .lphy files
        base_name=$(basename "$lphy_file" .lphy)
        slurm_script="${dir}/run_${base_name}.sh"
        
        # Create SLURM job script
        cat > "$slurm_script" <<EOL
#!/bin/bash
#SBATCH -J SM_${base_name}
#SBATCH -A uoa04397
#SBATCH --time=2:00:00
#SBATCH --mem=2G
#SBATCH -o ${dir}/${base_name}.out
#SBATCH -e ${dir}/${base_name}.err

${lphybeastPath}/lphybeast -l 40000000 "${lphy_file}"
EOL

        echo "Created SLURM script: $slurm_script"
        #sbatch "$slurm_script"
    done
done

echo "All SLURM job scripts generated."