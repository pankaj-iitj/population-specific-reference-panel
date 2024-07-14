## bash script for genotype quality control
## input files expected to be in VCF/ directory and named like chr<n>.vcf
## output will plink binary files (bed,bim,fam), pedigree (ped) and map (.map) files

for n in {1..22}; do 
  mkdir -p plink_binary/chr$n 
  plink2 --vcf vcf/chr$n.vcf --double-id --export ped --make-bed --out plink_binary/chr$n/chr$n 
done
