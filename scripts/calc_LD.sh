#!/bin/bash

## first step is to extract blocks from the ranges in bed files
# Loop over chromosome numbers from 1 to 22
for i in {1..22}
do
  BED_FILE="ld_blk/chr${i}.bed"
  PLINK_BFILE="geno_qc/chr${i}"
  OUTPUT_DIR="ld_ref/chr${i}/blocks"
  mkdir -p $OUTPUT_DIR
  
  # Initialize the block counter
  num=1
  
  # Loop over each line in the BED file
  while IFS=$'\t' read -r chromosome start end
  do
    # Skip the header line
    if [ "$chromosome" != "chromosome" ]; then
      # Execute the PLINK command to extract the specified block region
      plink --bfile $PLINK_BFILE --chr $i --from-bp $start --to-bp $end --make-bed --out $OUTPUT_DIR/block${num}
      # Increment the block number
      ((num++))
    fi
  done < $BED_FILE

  echo "Extraction completed for chromosome $i!"
done

echo "All extractions completed!"

# Next step to calculate LD matrix for each block in each chromosome

echo "Calculating LD matrix" 
for chr in {1..22}
do
  # Get the number of blocks for the current chromosome by counting the .bim files in the directory
  num_blocks=$(ls ld_ref/chr${chr}/blocks/block*.bim | wc -l)
  
  # Loop through each block for the current chromosome
  for n in $(seq 1 $num_blocks)
  do
    # PLINK command for each block
    plink --keep-allele-order --bfile ld_ref/chr${chr}/blocks/block${n} --r --matrix --out led_ref/chr${chr}/blocks/block${n}
  done
done
