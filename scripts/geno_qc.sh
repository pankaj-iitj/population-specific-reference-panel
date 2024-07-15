## bash script for quality control of plink binary genotype data
## input expects to be in plink_binary/ 
## output clean binaries to be in geno_qc/ directory. Binaries will be segregated per chromosome

# Loop through each chromosome
for n in {1..22}; do
  mkdir -p geno_qc/chr$n
  plink2 --bfile plink_binary/chr$n/chr$n --mind 0.01 \
  --geno 0.01 \
  --maf 0.01 \
  --hwe 1e-6 \
  --het \
  --make-bed \
  --out geno_qc/chr$n/chr${n}.clean1
done

# Loop through each chromosome for excludes samples with heterozygosity (+-3SD i.e. F=+-0.1)
for n in {1..22}; do
  awk 'NR>1 && $6>0.1 || $6<-0.1 {print $1,$2}' geno_qc/chr$n/chr${n}.clean1.het > geno_qc/chr$n/chr${n}_high_het.sample 
  plink2 --bfile geno_qc/chr$n/chr${n}.clean1 --remove geno_qc/chr$n/chr${n}_high_het.sample --make-bed --out geno_qc/chr$n/chr${n}.clean2
  
  rm geno_qc/chr$n/chr${n}.clean1.het geno_qc/chr$n/chr${n}.clean2.log
  mv geno_qc/chr$n/chr${n}.clean2.bim geno_qc/chr$n/chr${n}.bim
  mv geno_qc/chr$n/chr${n}.clean2.bed geno_qc/chr$n/chr${n}.bed
  mv geno_qc/chr$n/chr${n}.clean2.fam geno_qc/chr$n/chr${n}.fam
done
