#!/bin/bash
#SBATCH -J SM_runLphybeast
#SBATCH -A nesi00390
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --hint=multithread
#SBATCH --mail-user=yil218@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
# Define lphybeast executable and version files
lphybeastPath="/home/yil218/beast/bin"
phyloncoVersionFile="/home/yil218/.beast/2.7/phylonco/version.xml"
lphybeastVersionFile="/home/yil218/.beast/2.7/lphybeast/version.xml"
flcVersionFile="/home/yil218/.beast/2.7/lphybeast.flc/version.xml"

# Define your target directories
directories=( "./ErrModelScripts" "./withoutErrModelScripts" "./TrueErrModelScripts" "./TrueWithoutErrModelScripts" )

# Loop through each directory
for dir in "${directories[@]}"; do
    for lphy_file in "${dir}"/*.lphy; do
        [ -e "$lphy_file" ] || continue  # skip if no .lphy files
        base_name=$(basename "$lphy_file" .lphy)
        slurm_script="${dir}/run_${base_name}.slurm"
        
        # Create SLURM job script
        cat > "$slurm_script" <<EOL
#!/bin/bash
#SBATCH -J SM_${base_name}
#SBATCH -A nesi00390
#SBATCH --time=16:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --hint=multithread
#SBATCH --mail-user=yil218@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH -o ${dir}/${base_name}.out
#SBATCH -e ${dir}/${base_name}.err

module load beast

${lphybeastPath}/lphybeast -vf "${phyloncoVersionFile},${lphybeastVersionFile},${flcVersionFile}" -l 20000000 "${lphy_file}" > "${lphy_file}.log" 2>&1
EOL

        echo "Created SLURM script: $slurm_script"
        sbatch "$slurm_script"
    done
done

echo "All SLURM job scripts generated."
