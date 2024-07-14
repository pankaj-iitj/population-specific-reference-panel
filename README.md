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

### Joint calling:
After all the trimmed reads are deposited in the `trimmed_reads` directory, run `bash joint_calling.sh` command to execute joint calling step over multiple samples. The script takes trimmed reads from each sample and then aligns them with the reference genome, marks duplicate reads, and recalibrates base calls and calls variants per sample before finally perform joint calling over all samples. These include the commands: `bwa mem`, `samtools view`, `samtools sort`, `samtools index`, `gatk4 MarkDuplicatesSpark`, `gatk4 BaseRecalibrator`, `gatk4 ApplyBQSR`, `gatk4 HaplotypeCaller`, `gatk4 GenomicsDBImport` and `gatk4 GenotypeGVCFs` one after another. 

**Note:** We parallelized this process with the GNU `parallel` command taking 10 concurrent jobs and adjusted the number of threads as per our in-house server. Users are advised to calculate these parameters according to the memory allocation their system permits. 

The joint calling will yield multi-sample VCF files in `vcf/` directory split over chromosomes. During the `GenotypeGVCFs` execution, we only included 22 autosomes while sex-chromosomes and mDNA were excluded. ***This step is computationally intensive and may take several days depending on file sizes and number of samples***. 

### LD-block partitioning: 
The first step here is to convert the genotyped VCF files into PLINK binary (bed,bim,fam) files. The pipeline expects PLINK files to be sub-directories specified by each chromosome inside the `ld_blk/` directory. To do this run the command below
```
for n in {1..22}; do 
  mkdir -p ld_blk/chr$n 
  plink2 --vcf vcf/chr$n.vcf --double-id --export ped --make-bed --out ld_blk/chr$n/chr$n 
done
```
***Caution:*** It is quite common that during file conversion to plink binaries, allele swap happens. PLINK2 treats **minor allele as A1 and major allele as A2**. The current release of PLINK2 does not affect on `--keep-allele-order` flag as it handles A1/A2 alleles separately from **'reference allele'**. If you are still using PLINK1.9, use ``plink --keep-allele-order`` option to prevent allele swap. To check whether there had been allele swap or strand flipping we recommend using the commands beforehand:
```
bcftools +fixref chr<n>.vcf -- -f $genome # check for REF/ALT flip
snpflip --fasta-genome=$genome --bim-file=chr<n> -o chr<n> # correct the swap if needed
```
In the next step, we provided an R-script `interpolate.R` to interpolate the recombination rates (in cM) from known genetic distances between HapMap variants. This is required as during the conversion of plink binaries from VCF the genetic distance (cM) information in the *.map* file is lost. We used known recombination rates for the 1000G **South Asian (SAS)** populations available at [Pyrho recombination map](https://github.com/popgenmethods/pyrho?tab=readme-ov-file#human-recombination-maps) repository. Run `R --vanilla > interpolate.R` command to interpolate the genetic distances between the variants from VCF files to tab-separated text files. The first column of the output file is chromosome name, second column is the position of the variant in bp and the third column is the recombination rate in cM which look like this:
```
chr1	16103	0.000522530060230873
chr1	51479	0.00167045426135659
chr1	51898	0.00168405049157684
chr1	51928	0.00168502396868092
chr1	51954	0.00168586764883778
chr1	54490	0.00176815891336896
chr1	54669	0.0017739673267566
chr1	54708	0.0017752328469919
```










