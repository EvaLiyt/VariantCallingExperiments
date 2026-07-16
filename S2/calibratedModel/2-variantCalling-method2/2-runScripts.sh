for file in scripts/run_r*sh;
do
sbatch ${file}
done
