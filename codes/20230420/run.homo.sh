#!/bin/bash -l

# check how many homozygous of rs6910507 in our GWAS samples
module load htslib/1.16 bcftools/1.16

# batch1
batch1_sample_dir='/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/GWAS_UNITE/endophenotypes/batch1/CA1/chr6/sample.txt'
batch1_gwas_dir='/restricted/projectnb/cte/Challenge_Project_2022/data/genetic/1st_batch/_TOPMed_imputation/chr6.dose.vcf.gz'
bcftools query -r chr6:29995149 -f '%CHROM\t%POS[\t%GT]\n' -S ${batch1_sample_dir} ${batch1_gwas_dir} > rs6910507.homo.txt

# batch2
batch2_sample_dir='/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/GWAS_UNITE/endophenotypes/batch2/CA1/chr6/sample.txt'
batch2_gwas_dir='/restricted/projectnb/cte/Challenge_Project_2022/data/genetic/2nd_batch/_TOPMed_imputation/chr6.dose.vcf.gz'
bcftools query -r chr6:29995149 -f '%CHROM\t%POS[\t%GT]\n' -S ${batch2_sample_dir} ${batch2_gwas_dir} >> rs6910507.homo.txt

count=`grep -o '1|1' rs6910507.homo.txt | wc -l`
echo "${count} of them are homozygotes with minor allele G"
