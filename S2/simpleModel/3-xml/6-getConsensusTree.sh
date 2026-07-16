#!/bin/bash

#SBATCH -J CM_getConsensusTrees
#SBATCH --time=12:00:00        # Adjust walltime as needed


# Run the Python script
cd $HOME/beast-2.7/bin
folder="${WORK_DIR:?set WORK_DIR to your own working directory path}/vcExp/simpleModel/3-xml/results"

for filename in ${folder}/*trees; do
    echo "processing ${filename}"
  ./treeannotator -topology CCD0 -height CA ${filename} "${filename}.txt"
  echo "done ${filename}"
done
