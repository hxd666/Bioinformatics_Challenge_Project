#extract SNPs from chr17:45387303 -- 47117396
module load htslib/1.16 bcftools/1.16

printf "CHR\tBP\tALT\tREF\n" > TOPMed.CTE.hg38.snps
bcftools query -r chr17:45387303-47117396 -f '%CHROM\t%POS\t%ALT\t%REF\n' /restricted/projectnb/cte/Challenge_Project_2022/data/genetic/1st_batch/_TOPMed_imputation/chr17.dose.vcf.gz >> TOPMed.CTE.hg38.snps
