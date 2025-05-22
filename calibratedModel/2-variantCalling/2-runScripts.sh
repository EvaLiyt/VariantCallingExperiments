#!/bin/bash
#SBATCH -J CM_submitJobs
#SBATCH -A nesi00390            # The account code
#SBATCH --time=2:00:00        # Adjust walltime as needed
#SBATCH --mem=2G              # Adjust memory as needed
#SBATCH --cpus-per-task=4      # Number of CPUs per task
#SBATCH --mail-user=yil218@aucklanduni.ac.nz
#SBATCH --mail-type=ALL

for i in {0..99}; do
    repeat="r$i"
    script_path="../1-simulation/pileup/${repeat}/run_${repeat}.sh"
    
    if [[ -f "$script_path" ]]; then
        echo "Submitting $script_path"
        chmod +x "$script_path"
        sbatch "$script_path"
    else
        echo "Script not found: $script_path"
    fi
done
