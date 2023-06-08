module load htslib/1.16 bcftools/1.16

printf "CHR\tBP\tALT\tREF\n" > GRCh38_alnloci_mpileup_chrm.snps
bcftools query -f '%CHROM\t%POS\t%ALT\t%REF\n' /restricted/projectnb/cte/Challenge_Project_2022/jnpetros/dbSNP/redo/GRCh38_alnloci_mpileup_chrm.vcf >> GRCh38_alnloci_mpileup_chrm.snps
