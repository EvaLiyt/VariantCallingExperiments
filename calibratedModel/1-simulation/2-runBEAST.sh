module load Python/3.11.3-gimkl-2022a
python modifyXML.py

file="./modified.xml"
bash temp-beast.sh ${file}
sbatch "${file}.sh"
