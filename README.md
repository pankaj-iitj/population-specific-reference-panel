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
| [BWA](https://bio-bwa.sourceforge.net/) | Read mapping to reference genome |
| [Samtools](https://www.htslib.org/) | Sorting and indexing alignment files |
| [GATK v4.6.0](https://gatk.broadinstitute.org/hc/en-us) | Toolkit for various functions e.g. mark duplicates, base recalibration, variant calling |
| [PLINK v2.0](https://www.cog-genomics.org/plink/2.0/) | Genotype file conversion, LD calculation |
| [bcftools](https://samtools.github.io/bcftools/howtos/index.html) | VCF file manipulation |      
| [snpflip](https://github.com/biocore-ntnu/snpflip) | Correcting strand flipping and allele swap |
| [LDetect](https://bitbucket.org/nygcresearch/ldetect/src/master/) | LD Block partitioning | 

## Set the working directory
We assume per sample paired-end reads i.e. `*R1.fastq` and `*R2.fastq` are in the `reads/` directory. The required human genome reference (GRCh37/38) in FASTA and dnSNP [build 156](https://ftp.ncbi.nih.gov/snp/archive/b156/VCF/) known sites in VCF should be placed in the `reference/` directory. We provided the `download_ref.sh` file for manual download in case the user has a choice of using their preferred reference version. Please make directories for intermediate steps by the following command: 
```
mkdir quality trimmed_reads sam bam sorted_bam dupl_bam bqsr_bam gvcf vcf geno_qc ld_blk ld_ref
```
The `joint_call.sh` script also creates the directories, so if you decide simply to run the joint calling script, you can skip the above command. 

## Prepare input data
The only initial step needed is to index the reference genome by `bwa index` and `samtools faidx` as they are needed while running these commands during read mapping and alignment sorting. For joint calling, an additional step of creating a sequence dictionary is needed using `GATK CreateSequenceDictionary`. 
```
genome_dir="reference"
genome='reference/<ref>.fa    # point towards reference fasta file
bwa index $genome
samtools faidx $genome
gatk4 CreateSequenceDictionary R=$genome O=$genome_dir/<ref>.dict
```
## Run the scripts
### Quality check before joint calling:
Run the `quality.sh` script to perform the steps one by one: quality check by Fastqc, multi-sample aggregated quality report by MultiQC and read trimming by fastp. If read qualities are okay, then no need for trimming and in that case, comment out the fastp command line:
```
cat samples.txt | parallel --progress --eta -j 10 "fastp -i reads/{}_R1.fastq.gz -I reads/{}_R2.fastq.gz -o trimmed_reads/{}_R1.trim.fastq.gz -O trimmed_reads/{}_R2.trim.fastq.gz"
```
The quality reports will be dumped in the `quality` directory. 

### joint calling:










