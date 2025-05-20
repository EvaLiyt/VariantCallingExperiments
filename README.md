# VariantCallingExperiments
This is a variant calling simulation study pipeline, including realistic sequences simulation, single cell sequencing simulation, variant calling, and Bayesian analysis. For simplicity, please put all the software in the same folder "Workspace".

## Environment
- **Operating System**: MacOS Ventura v13.4.1
- **Java Version**:
  - OpenJDK version "17.0.12" (LTS), released on July 16, 2024
  - OpenJDK Runtime Environment Zulu 17.52+17-CA (build 17.0.12+7-LTS)
  - OpenJDK 64-Bit Server VM Zulu 17.52+17-CA (build 17.0.12+7-LTS, mixed mode, sharing)
- **Python 3**
- **Python 2.7**: only for SCcaller

## Materials
* [NCBI Datasets command-line tools](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/)
* [LinguaPhylo v1.6.1](https://linguaphylo.github.io/)
* [SCcaller v2.0.0](https://github.com/biosinodx/SCcaller)
* [BEAST2 v2.7.6](https://www.beast2.org/)

## Pipeline

### Clone This Repo
1. Open your terminal and navigate to your working directory.
2. Clone the repository(repo) to your local with the following code:
   ```
   git clone https://github.com/EvaLiyt/VariantCallingExperiments.git
   ```
### Preparation
1. Download human reference genome GRCh38.p14 from NCBI datasets.
2. Obtain 1 million sites from p-arm with `obtainSequence.py`.
3. Calculate parameters with `priorCalc.Rmd`. (TODO: wait to finish)
4. 
