#!/bin/bash

# bash script to automate partitioning chromosome into nearly independent LD-blocks
# assumes input clean variant data is in geno_qc/chr${n}/ directory
# e.g. for chromosome 22, input files will be geno_qc/chr22/chr22.bed, geno_qc/chr22/chr22.bim and geno_qc/chr22/chr22.fam
# assumes Ldetect already installed (see, 'install.sh')

# activate 'ldetect' virtual environment
source ~/ldetect/bin/activate
PYTHON_VERSION=$(basename $(ls -d ~/ldetect/lib/python* | head -n 1))
LDETECT_PACKAGES_DIR=~/ldetect/lib/$PYTHON_VERSION/site-packages/ldetect
TARGET_DIR=$PWD/ld_blk
INDIR=$PWD/geno_qc

if [ -d "$LDETECT_PACKAGES_DIR" ]; then
    echo "Copying ldetect module from $LDETECT_PACKAGES_DIR to $TARGET_DIR"
    cp -r "$LDETECT_PACKAGES_DIR" "$TARGET_DIR"
    echo "DONE"
else
    echo "Ldetect module not found! Abort! ! !"
fi
# point working directory
echo "Pointing towards base working directory for Ldetect"
BASEWD=$TARGET_DIR/ldetect/examples
cd $BASEWD
echo "We are now at base directory. We will work from here for a while"
# make output directories
mkdir example_data
mkdir example_data/cov_matrix
mkdir example_data/cov_matrix/scripts
mkdir example_data/cov_matrix/chr22
mkdir example_data/vector example_data/minima example_data/bed
# prepare input clean genotype files from 'geno_qc' dir
mkdir VCF
plink --bfile $INDIR/chr22/chr22 --keep-allele-order --recode vcf --out VCF/chr22
bgzip -c VCF/chr22.vcf > VCF/chr22.vcf.gz
tabix -p vcf VCF/chr22.vcf.gz
rm VCF/chr22.log VCF/chr22.nosex

echo "Input files prepared !"

##########################################################################################################################################################################################################
# prepare genotyped sample panel from 1000 Genomes Project Phase 3 call-set
echo "Preparing sample set .... "
bcftools query -l VCF/chr1.vcf.gz > example_data/sample.txt
N_SAMPLE=$(wc -l < example_data/sample.txt)
echo "Sample set ready, total sample: $N_SAMPLE"

#########################################################################################################################################################################################################
# prepare genetic map file
# make sure genetic map is in example_data/genetic_map/chr22.tab.gz
echo "Preparing genetic map file .... "
file_path="example_data/genetic_map/chr22.tab.gz"
if [[ ! -f $file_path ]]; then
   echo "Please move genetic map files into example_data/genetic_map/ directory with corresponding chromosomes."
   exit 1  # Terminate the program with an error code
fi

##########################################################################################################################################################################################################
# Calculate LD breakpoints
# Step 1: partition chromosome
genetic_map="example_data/genetic_map/chr22.tab.gz"
n_sample=417
output_file="example_data/cov_matrix/scripts/chr22_partitions" 
python3 P00_00_partition_chromosome.py $genetic_map $n_sample $output_file
echo "partition done"
echo "Computing covariance matrices"

# Step 2: compute covariance matrix
awk '{print "tabix -h VCF/chr22.vcf.gz 22:" $1 "-" $2 " | python3 P00_01_calc_covariance.py example_data/genetic_map/chr22.tab.gz example_data/sample.txt 11418 1e-7 example_data/cov_matrix/chr22/chr22." $1 "." $2 ".gz"}' chr22_partitions > commands.txt
parallel --progress -j 18 < commands.txt
rm commands.txt
echo "Covariance matrix calculated ! Proceeding to step next step"

# Step 3: convert covariance matrix to vector
python3 P01_matrix_to_vector_pipeline.py --dataset_path=example_data/cov_matrix/ --name=chr22 --out_fname=example_data/vector/vector-SAS-chr22.txt.gz
echo "Done ! Now converting into vector"

# Step 4: calculate minima
python3 P02_minima_pipeline.py \
        --input_fname=example_data/vector/vector-SAS-chr22.txt.gz \
        --chr_name=chr22 --dataset_path=example_data/cov_matrix/ \
        --n_snps_bw_bpoints=7000 --out_fname=example_data/minima/minima-SAS-chr22.pickle
echo "minima calculated successfully"

# Step 5: extract final breakpoints
python3 P03_extract_bpoints.py \
        --name=chr22 --dataset_path=example_data/cov_matrix/ \
        --subset=fourier_ls \
        --input_pickle_fname=example_data/minima/minima-SAS-chr22.pickle > example_data/bed/BED-SAS-chr22.bed
echo "Uniform breakpoints are at example_data/bed/ directory"



