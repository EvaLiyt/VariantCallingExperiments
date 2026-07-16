#!/bin/bash
#SBATCH -J RM_submitJobs
#SBATCH --time=2:00:00        # Adjust walltime as needed

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
