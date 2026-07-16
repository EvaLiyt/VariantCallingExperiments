# VariantCallingExperiments
This is a variant calling simulation study pipeline, including realistic sequences simulation, single cell sequencing simulation, variant calling, and Bayesian analysis. For simplicity, please put all the software in the same folder "Workspace".

## Environment
- **Operating System**: MacOS Ventura v13+
- **Java Version**:
  - OpenJDK version "17.0.12" (LTS), released on July 16, 2024
  - OpenJDK Runtime Environment Zulu 17.52+17-CA (build 17.0.12+7-LTS)
  - OpenJDK 64-Bit Server VM Zulu 17.52+17-CA (build 17.0.12+7-LTS, mixed mode, sharing)
- **Python 3**

## Materials
* [NCBI Datasets command-line tools](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/)
* [LinguaPhylo v1.7.0](https://linguaphylo.github.io/)
* [SCcaller v2.0.0](https://github.com/biosinodx/SCcaller)
* [BEAST2 v2.7.8](https://www.beast2.org/)

## SLURM account and personal paths
All SLURM job scripts in this pipeline read the submitting account from the standard
`SBATCH_ACCOUNT` environment variable rather than hardcoding one, so before submitting
any job set your own account once per session (or add it to your shell profile):
```
export SBATCH_ACCOUNT=your_account_code
```
A few scripts also build absolute paths into your own storage (results, merged VCFs). These
read a separate `WORK_DIR` variable holding your full personal base path, wherever that is
on your system:
```
export WORK_DIR=/path/to/your/own/working/directory
```
Both variables will fail with a clear error message if left unset when a script needs them.
No script needs to be edited. Paths under your own `$HOME` (installed software, reference
data) are resolved automatically and need no configuration.

## Preparation
### Clone The Repos
1. Open your terminal and navigate to your working directory.
2. Clone the repository(repo) to your local with the following code:
   ```
   git clone https://github.com/EvaLiyt/VariantCallingExperiments.git
   ```

### Generate Lphy Scripts
1. Download human reference genome GRCh38.p14 from NCBI datasets.
	```
	datasets download genome accession GCF_000001405.40 --filename human_GRCh38_dataset.zip
	unzip human_GRCh38_dataset.zip
	```
	This unzips into `ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna`.
	Create `~/data/` if it does not already exist, then move that `.fna` file there so the path
	`~/data/GCF_000001405.40_GRCh38.p14_genomic.fna` exists exactly as named. `1-obtainSeq.py`
	reads this exact path.
2. Obtain a 10,000 bp region from the p-arm with `1-obtainSeq.py`. This writes `~/data/10K.fna`,
	which every downstream script (in every language) reads from that same path.
3. Get the homozygous reference sequence as well with `2-getHomoRoot.py`. This writes
	`~/data/10K_homo.fna`.
4. Calculate parameters and get Lphy scripts with `S1/getLphyScripts.R`. (Or the 40 tumour cells
	simulation with `S2/getMoreTaxaScripts.R`)

## Read Count Simulation
### Simple Model
Simple model is a 40 tumour cells Yule tree with same branch rate.
1. Run simulation for a hundred times with `1-runSimulation.sh`. The results will be annotated and allocate to different folders.
2. Convert the read count data to pileup files with `2-runRC.R`.

### Realistic Model
Realistic model is a subsampled tree from 2000 taxa trees containing normal and tumour cells. The branch rates are different between normal and tumour cells.
1. Run simulation for a hundred times with `1-runSimulation.sh`. The results will be annotated and allocate to different folders.
2. Convert the read count data to pileup files with `2-runRC.R`.

### Calibrated Model
This will be finished when we have a right simulator (will use calibratedCPP then).

### Variant Calling 
SCCaller is modified to accept pileup files since read count data does not have features to create bam file. Python 2.7 is also upgraded to Python 3+ for future use.
Vcf2fasta is modified in the fork, will be merge to the main branch when it is satisfied by the author.
1. Generate a random set of minDepth and minvar for each dataset to fill the slurm job script with `1-generageScript.py`. Each set of parameters are written down in a separate file. Run them with `2-runScripts`.
2. Filter and adjust vcf formats to merge each dataset with `3-runFollowing`.
3. Use `vcf2fasta` to convert merged vcf to GT16 alignment with `4-generateMappedFile.sh`.
4. Replace uncertainties in the final alignment with homozygous reference sequence with `5-replaceQuestionMarks.sh`.

### Bayesian Inference
1. Get slurm jobs for each lphyBeast script with '1-getSlurmJobs.sh'.
2. Check all the xmls are generated, then use '2-modifyXMLs.py' to replace the operators from 'phylonco' to suit our data.
3. Run '3-getSH.sh' to get bash scripts then use '4-runBEAST.sh' to submit BEAST2 runs.
4. Organise the result files with '5-organiseResults.sh'.

### MethodComparison
The `methodComparison` folder contains methods comparison scripts for different ways to map variant called results.

Method1 is `2-variantCalling`, which uses the strong evidence and map missing sites with root.
Method2 is `2-variantCalling-method2`. In here we use the D1V1 files to map the sites in the merged VCF but with no genotype. Other missing sites mapped with root.
Method3 is `2-variantCalling-method3`. In here we use the D1V1 files to map the sites in the merged VCF but with no genotype.
Method4 is `2-variantCalling-method4`. Variant sites are filtered with minD2 minVar0 and pass back in `sccaller_labVersion_2ndRound.py` for variant calling.



