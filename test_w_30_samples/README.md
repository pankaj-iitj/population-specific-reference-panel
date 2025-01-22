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

Make sure the files are correctly placed in each directory. 
