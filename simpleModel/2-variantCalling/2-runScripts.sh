#!/bin/bash
#SBATCH -J SM_submitJobs
#SBATCH -A uoa04397 
#SBATCH --time=2:00:00 

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
