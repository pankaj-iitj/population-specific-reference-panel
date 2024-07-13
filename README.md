# Comprehensive pipeline for creating LD-based reference panel
![Pipeline_workflow](https://github.com/user-attachments/assets/a8143c56-f936-4ceb-8f7f-c67a80478613)

We present a comprehensive pipeline to creating Linkage Disequilibrium (LD)-based reference panel from whole genome sequencing (WGS) paired-end reads. The pipeline expects whole genome paired-end reads (R1/R2) in `.fastq (or, .fastq.gz)` and the final output will be an LD-based reference panel in HDF5 format. 
Broadly, the pipeline is divided into three main phases:
* Joint calling phase
* LD-block partitioning phase
* LD calculation and reference panel creation phase

Below, we explain how to start with the required tools, set the working directory, prepare input files and run the scripts.
## Required tools
The tools included in this pipeline are listed sequentially. Make sure these tools are pre-installed. If not, please run the `install_tools.sh` file to install them beforehand. 

| Tools | Description |
| --- | --- |
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | Quality control of paired-end `fastq' files |
| [MultiQC](https://multiqc.info/) | Aggregated quality report |
| [fastp](https://github.com/OpenGene/fastp) | Read trimming in batch |
