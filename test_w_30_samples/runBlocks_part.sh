#!/bin/bash
mkdir example_data && cd example_data
mkdir bed cov_matrix minima vector && cd cov_matrix
mkdir scripts
for i in {1..22}; do mkdir chr${i}; done && cd ../../
# Step: 1 (partition chromosome)
for i in {1..22}; do python3 P00_00_partition_chromosome.py example_data/chr${i}.genetic_map.tab.gz 30 example_data/cov_matrix/scripts/chr${i}_partitions; done
# Step: 2 (calculate covariance)
for chr in {1..22}; do
    partition_file="example_data/cov_matrix/scripts/chr${chr}_partitions"
    if [[ -f $partition_file ]]; then
        awk -v chr=$chr '{print "tabix -h qc_vcf/chr" chr ".vcf.gz chr" chr ":" $1 "-" $2 " | python3 P00_01_calc_covariance.py example_data/chr" chr ".genetic_map.tab.gz example_data/sample.txt 14296 1e-7 example_data/cov_matrix/chr" chr "/chr" chr "." $1 "." $2 ".gz"}' $partition_file > commands_chr${chr}.txt
    else
        echo "Partition file not found for chr${chr}: $partition_file"
    fi
done
for i in {1..22}.txt; do
    parallel --progress -j 20 < commands_chr${i}.txt
done
# Step: 3 (convert covariance matrix to vector)
for i in {1..22}; do
    python3 P01_matrix_to_vector_pipeline.py \
    --dataset_path=example_data/cov_matrix/ --name=chr${i} --out_fname=example_data/vector/vector.chr${i}.txt.gz
done
# Step: 4 (calculate minima)
for i in {1..22}; do
    pythob3 P02_minima_pipeline.py \
    --input_fname=example_data/vector/vector.chr${i}.txt.gz --chr_name=chr${i} \
    --dataset_path=example_data/cov_matrix/ \
    --n_snps_bw_bpoints=7000 --out_fname=example_data/minima/minima.chr${i}.pickle
done
# Step: 5 (extract breakpoints)
for i in {1..22}; do 
    python3 P03_extract_bpoints.py --name=chr${i} --subset=fourier_ls \
    --dataset_path=example_data/cov_matrix/ \ 
    --input_pickle_fname=example_data/minima/minima.chr${i}.pickle > example_data/bed/chr${i}.bed 
done
