#!/bin/bash

# Initialize reference directories
genome_dir="reference/"
genome="reference/hg38.fa"
# head -n 10 $genome

known_sites_dir="reference/"
known_sites="reference/dbsnp151-00-All.vcf.gz"
# zcat $known_sites | head

reads="reads"         

# bwa index $genome                                          [                                ]
# (or),                                                      [   Indexing Reference Genome    ]  
# samtools faidx $genome                                     [                                ]


# gatk CreateSequenceDictionary R=$genome O=$genome_dir/hg38.dict                              [ Create Dictionary for Reference Genome ] 


############################################################################### < Main Workflow > #########################################################################################################


# mkdir sam bam sorted_bam dupl_bam bqsr_bam gvcf vcf

for fq1 in $reads/*_R1_001.fastq;
do
    echo "working with file $fq1"; base=$(basename $fq1 _R1_001.fastq);  echo "base name is $base" \

    fq1=$reads/${base}_R1_001.fastq
    fq2=$reads/${base}_R2_001.fastq
   
    sam=sam/${base}.aligned.sam
    bam=bam/${base}.aligned.bam
    
    sorted_bam=sorted_bam/${base}.aligned.sorted.bam
    dupl_bam=dupl_bam/${base}.aligned.sorted.dupl.bam
    bqsr_bam=bqsr_bam/${base}.aligned.sorted.bqsr.bam
    gvcf=gvcf/${base}_variants.g.vcf
    final_variants=vcf/${base}_final_variants.vcf

    bwa mem -t 20 -R "@RG\tID:$base\tPL:ILLUMINA\tSM:$base" $genome $fq1 $fq2 > $sam
    samtools view -@ 20 -S -b $sam > $bam
    samtools sort -@ 20 -o $sorted_bam $bam
    samtools index $sorted_bam
    gatk MarkDuplicatesSpark -I $sorted_bam -O $dupl_bam -M dupl_bam/${base}_Mark_dupl_Metrics.txt --conf 'spark.executor.cores=16'
    gatk BaseRecalibrator -I $dupl_bam -R $genome --known-sites $known_sites -O dupl_bam/${base}_recall_data.table
    gatk ApplyBQSR -I $dupl_bam -R $genome --bqsr-recal-file dupl_bam/${base}_recall_data.table -O $bqsr_bam
    gatk --java-options "-Xmx4g" HaplotypeCaller -R $genome -I $bqsr_bam -O $gvcf -ERC GVCF --dbsnp $known_sites
done    
# List of chromosomes (adjust if needed)
chromosomes=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
             chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 \
             chr21 chr22 chrX chrY chrM)

# Loop over chromosomes
for chr in "${chromosomes[@]}"; do
    echo "Processing $chr"
    
    # GenomicsDBImport
    gatk GenomicsDBImport --java-option "Xmx2G Xms2G" -R $genome --genomicsdb-workspace-path "vcf/gendb/${chr}"  -L "$chr"  $(for gvcf in gvcf/*.g.vcf; do echo -n "-V $gvcf "; done) --batch-size 50 --reader-threads 5

    # GenotypeGVCFs
    gatk GenotypeGVCFs --java-option "Xmx2G Xms2G" -R "$genome" --gendb gendb/${chr} -O "vcf/${chr}_final.vcf"
done
