#!/bin/sh

echo "#!/bin/sh" > $1.sh
echo "" >> $1.sh
echo "#SBATCH -J $1		# The job name" >> $1.sh
echo "#SBATCH --time=2:00:00         # The walltime" >> $1.sh
echo "#SBATCH --mem=4G 	        # in total" >> $1.sh
echo "#SBATCH --cpus-per-task=2       # OpenMP Threads" >> $1.sh
echo "#SBATCH --ntasks=1              # not use MPI" >> $1.sh
echo "#SBATCH --hint=multithread      # A multithreaded job, also a Shared-Memory Processing (SMP) job" >> $1.sh
echo "#SBATCH -D ./			# The initial directory" >> $1.sh
echo "#SBATCH -o $1_out.txt		# The output file" >> $1.sh
echo "#SBATCH -e $1_err.txt		# The error file" >> $1.sh

echo "# sacct -j JOBID --format=\"ReqMem,MaxRSS,CPUTime,AveCPU,Elapsed\"" >> $1.sh

echo "module load beagle-lib/4.0.0-GCC-11.3.0" >> $1.sh
echo "module load Java/17" >> $1.sh

echo "# beast 2.7.x " >> $1.sh
echo "srun $HOME/beast-2.7/bin/beast -beagle_SSE $1" >> $1.sh
