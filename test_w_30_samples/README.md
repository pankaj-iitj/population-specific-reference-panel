Test block partitioning with a small sample (min. 30) size. Here, we assume you have a QCed genotype in the VCF file. For convenience, we gave an idea of the directory hierarchy as it should be during the computation. **Note:** For 22 autosomes, we recommend an HPC with at least 30 CPUs. Although it will work on a normal user laptop (e.g., 8 CPUs), the computation time will be significantly higher (~ 8 hrs).

# General instructions and file preparation
**1.** Follow the general directory convention for smooth running

**$PWD**
  * P00_00_partition_chromosome.py
  * P00_01_calc_covariance.py
  * P01_matrix_to_vector_pipeline.py
  * P02_minima_pipeline.py
  * P03_extract_bpoints.py
  * **example_data**
      * samples.txt
      * chr1.genetic_map.tab.gz
        
        .
        .
        
      * chr22.genetic_map.tab.gz
      * **cov_matrix**
          * **scripts**
          * **chr1**
            
            .
            .
            
          * **chr22**
      * **vector**
      * **minima**
      * **bed**
  * **qc_vcf**
      * chr1.vcf.gz
      * chr1.vcf.gz.tbi
        
        .
        .
        
      * chr22.vcf.gz
      * chr22.vcf.gz.tbi

Make sure the files are correctly placed in each directory. Directories are indicated with bold text.

**2.** Use `bcftools query -l *file.vcf.gz* > sample.txt` to list the sample ids. 

**3.** Activate the `ldetect` virtual environment by `source $PATH_to_ldetect/bin/activate` command before start computing.

**4.** You are ready to go and the final output will be `bed` files for respective chromosomes at the `bed/` directory. While using the bash script provided cross-check with file names. We left it as we used hoping the user would change them accordingly. We provided the BED file with high-LD blocks for 22 autosomes generated for the 1KG SAS population (30 samples randomly) using the genotype reference panel from [GRCh38-based 2018 biallelic SNV](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/)  

**5.** To make the complete LD panel refer to the main repo `README.md` for further steps
