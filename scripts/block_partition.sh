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
for n in {1..22}; do mkdir example_data/cov_matrix/chr${n}; done
mkdir example_data/vector example_data/minima example_data/bed
# prepare input clean genotype files from 'geno_qc' dir
# Loop through each chromosome
mkdir VCF
for n in {1..22}; do
    plink --bfile $INDIR/chr${n}/chr${n} --keep-allele-order --recode vcf --out VCF/chr${n}
    bgzip -c VCF/chr${n}.vcf > VCF/chr${n}.vcf.gz
    tabix -p vcf VCF/chr${n}.vcf.gz
    rm VCF/chr${n}.log VCF/chr${n}.nosex
done
echo "Input files prepared !"

##########################################################################################################################################################################################################
# prepare genotyped sample panel from 1000 Genomes Project Phase 3 call-set
echo "Preparing sample set .... "
# bcftools query -l VCF/chr1.vcf.gz > example_data/sample.txt
# $N_SAMPLE=$(wc -l < example_data/sample.txt)
echo "Sample set ready, total sample: $N_SAMPLE"

#########################################################################################################################################################################################################
# prepare genetic map file
echo "Preparing genetic map file .... "
for n in {1..22}; do
    file_path="example_data/genetic_map/chr${n}.tab.gz"
    if [[ ! -f $file_path ]]; then
        echo "Please move genetic map files into example_data/genetic_map/ directory with corresponding chromosomes."
        exit 1  # Terminate the program with an error code
    fi
done

########################################################################################################################################################################################################
# STEP 1: Partition the chromosome into overlapping bins based on genetic distances
echo "STEP 1: Creating partitions for chromosomes"
sasinds_file="example_data/sasinds.txt"
n_sample=$(wc -l < $sasinds_file)
for n in {1..22}; do
    genetic_map_file="example_data/genetic_map/chr${n}.tab.gz"
    output_file="example_data/cov_matrix/scripts/chr${n}_partitions"
    python3 P00_00_partition_chromosome.py $genetic_map_file $n_sample $output_file
done
echo "Partitioning DONE ! all 22 chromosomes"

#######################################################################################################################################################################################################
# STEP 2: Compute covariance matrix for each partition
for chr in {1..22}; do
    partition_file="example_data/cov_matrix/scripts/chr${chr}_partitions"
    if [[ -f $partition_file ]]; then
        awk -v chr=$chr '{print "tabix -h VCF/chr" chr ".vcf.gz " chr ":" $1 "-" $2 " | python3 P00_01_calc_covariance.py example_data/genetic_map/chr" chr ".tab.gz example_data/sasinds.txt 11418 1e-7 example_data/cov_matrix/chr" chr "/chr" chr "." $1 "." $2 ".gz"}' $partition_file > commands_chr${chr}.txt
        parallel --progress -j 60 < commands_chr${chr}.txt

        rm -f commands_chr${chr}.txt
    else
        echo "Partition file for chromosome ${chr} not found: ${partition_file}"
    fi
done
echo "Covariance matrix calculated ! Proceeding to step next step"

########################################################################################################################################################################################################
# STEP 3: Covert covariance matrices to vectors
# Loop through chromosomes 1 to 22
for chr in {1..22}; do
    partition_file="example_data/cov_matrix/scripts/chr${chr}_partitions"

    if [[ -f $partition_file ]]; then
        first_val=$(awk 'NR==1 {print $1}' $partition_file)
        last_val=$(awk 'END {print $2}' $partition_file)

        command="python3 P01_matrix_to_vector_pipeline.py --dataset_path=example_data/cov_matrix/ --name=chr${chr} --out_fname=example_data/vector/vector-SAS-chr${chr}-${first_val}-${last_val}.txt.gz"
        commands+=("$command")
    else
        echo "Partition file for chromosome ${chr} not found: ${partition_file}"
    fi
done

printf "%s\n" "${commands[@]}" | parallel --progress -j 22

#######################################################################################################################################################################################################
# STEP 4: Calculate minima from vectors
# Loop through chromosome 1 to 22
commands=()
for chr in {1..22}; do
    partition_file="example_data/cov_matrix/scripts/chr${chr}_partitions"
    if [[ -f $partition_file ]]; then
        first_val=$(awk 'NR==1 {print $1}' $partition_file)
        last_val=$(awk 'END {print $2}' $partition_file)
        input_file="example_data/vector/vector-EUR-chr${chr}-${first_val}-${last_val}.txt.gz"
        output_file="example_data/minima/minima-SAS-chr${chr}-7000-${first_val}-${last_val}.pickle"
        # construct command with number of SNPs between breakpoints 7000 (balance between hg19 and hg38 based LD panels)
        command="python3 P02_minima_pipeline.py --input_fname=${input_file} --chr_name=chr${chr} --dataset_path=example_data/cov_matrix/ --n_snps_bw_bpoints=7000 --out_fname=${output_file}"
        # add the command to the array
        commands+=("$command")
    else
        echo "Partition file for chromosome ${chr} not found: ${partition_file}"
    fi
done
printf "%s\n" "${commands[@]}" | parallel --progress -j 22

######################################################################################################################################################################################################
# STEP 5: Extract final breakpoints
# Final loop
commands_bpoint=()
for n in {1..22}; do
    partition_file="example_data/cov_matrix/scripts/chr${n}_partitions"
    if [[ -f $partition_file ]]; then
        first_val=$(awk 'NR==1 {print $1}' $partition_file)
        second_val=$(awk 'END {print $2}' $partition_file)
        input_pickle="example_data/minima/minima-SAS-chr${n}-7000-${first_val}-${second_val}.pickle"
        output_bed="example_data/bed/SAS-chr${n}-7000-${first_val}-${second_val}.bed"
        command_bpoint="python3 P03_extract_bpoints.py --name=chr${n} --dataset_path=example_data/cov_matrix/ --subset=fourier_ls --input_pickle_fname=${input_pickle} > ${output_bed}"
        commands_bpoint+=("$command_bpoint")
    else
        echo "Partition file for chromosome ${n} not found: ${partition_file}"
    fi
done
printf "%s\n" "${commands[@]}" | parallel -j 22



    

