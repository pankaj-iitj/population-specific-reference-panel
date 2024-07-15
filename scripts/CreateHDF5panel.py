import os
import h5py
import numpy as np
import re

# Function to read LD matrix and SNP list from files in a block folder
def read_block_data(ld_matrix_file, snplist_file):
    # Read LD matrix
    ld_matrix = np.loadtxt(ld_matrix_file)

    # Read SNP list
    with open(snplist_file, 'r') as f:
        snplist = [line.strip() for line in f]

    return ld_matrix, snplist

# Main function to create HDF5 file for a given chromosome
def create_hdf5(chromosome):
    blocks_dir = f'chr{chromosome}/blocks/'
    output_hdf5_file = f'ldblk_igv_chr{chromosome}.hdf5'

    with h5py.File(output_hdf5_file, 'w') as hdf5_file:
        for block_folder in sorted(os.listdir(blocks_dir)):
            if block_folder.startswith('block') and block_folder.endswith('.ld'):
                # Extract block number from the folder name using regular expression
                match = re.search(r'\d+', block_folder)
                block_number = match.group(0) if match else 'unknown'
                block_name = f'blk_{block_number}'

                # Check if the group already exists
                if block_name in hdf5_file:
                    print(f"Group {block_name} already exists. Skipping.")
                    continue

                # Construct file paths for LD matrix and SNP list
                ld_matrix_file = os.path.join(blocks_dir, f'block{block_number}.ld')
                snplist_file = os.path.join(blocks_dir, f'snplist.blk{block_number}.txt')

                # Ensure the files exist before attempting to read them
                if not os.path.isfile(ld_matrix_file) or not os.path.isfile(snplist_file):
                    print(f"Files for block {block_number} not found. Skipping.")
                    continue

                # Read block data
                ld_matrix, snplist = read_block_data(ld_matrix_file, snplist_file)

                # Create the group for the block
                block_group = hdf5_file.create_group(block_name)

                # Create subgroups and datasets within the block group
                block_group.create_dataset('ldblk', data=ld_matrix)
                block_group.create_dataset('snplist', data=np.array(snplist, dtype='S'))

# Loop through chromosomes and create HDF5 files for each
for chr_num in range(2, 22):
    create_hdf5(chr_num)
    print(f"Job done for chromosome {chr_num}!")

print("All jobs done !!!")



