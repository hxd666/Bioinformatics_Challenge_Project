#!/bin/bash -l

# This is the pipeline for the imputation of MAPT subhaplotypes for tomped vcf files using impute.v2 and shapeit

# The argument is 1st or 2nd
batch=${1}

echo "==============================================================================================================================================="
echo "Extract SNPs in MAPT region from vcf files  and Convert them to plink format"
echo "==============================================================================================================================================="
mkdir /restricted/projectnb/cte/Challenge_Project_2022/data/genetic/${batch}_batch/_mapt_region/
mkdir /restricted/projectnb/cte/Challenge_Project_2022/data/genetic/${batch}_batch/_mapt_region/bed_files
module load plink/2.00a2.3

plink2 --vcf /restricted/projectnb/cte/Challenge_Project_2022/data/genetic/${batch}_batch/_TOPMed_imputation/chr17.dose.vcf.gz --make-bed --extract /restricted/projectnb/cte/Challenge_Project_2022/data/genetic/reference/mapt/candidate.snps.hg38.txt --out /restricted/projectnb/cte/Challenge_Project_2022/data/genetic/${batch}_batch/_mapt_region/bed_files/MAPT.topmed

echo "==============================================================================================================================================="
echo "Converting the  GRCh38 ref in bim file to GRCh37 ref"
echo "==============================================================================================================================================="

module load R/3.6.2

Rscript --vanilla /restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20221027/convert.bim.R ${batch}

echo "==============================================================================================================================================="
echo " Prephasing using SHAPEIT"
echo "==============================================================================================================================================="

mkdir /restricted/projectnb/cte/Challenge_Project_2022/imputation/${batch}_batch/prephasing/

# check
/restricted/projectnb/cte/_xdhan/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -check \
        -B /restricted/projectnb/cte/Challenge_Project_2022/data/genetic/${batch}_batch/_mapt_region/bed_files/MAPT.topmed \
        -M /restricted/projectnb/cte/Challenge_Project_2022/data/genetic/reference/mapt/genetic_map_chr17_combined_b37.txt \
        --input-ref /restricted/projectnb/cte/Challenge_Project_2022/data/genetic/reference/mapt/_1000GP_Phase3/1000GP_Phase3_chr17.hap.gz /restricted/projectnb/cte/Challenge_Project_2022/data/genetic/reference/mapt/_1000GP_Phase3/1000GP_Phase3_chr17.legend /restricted/projectnb/cte/Challenge_Project_2022/data/genetic/reference/mapt/_1000GP_Phase3/1000GP_Phase3.sample \
        --output-log /restricted/projectnb/cte/Challenge_Project_2022/imputation/${batch}_batch/prephasing/gwas.alignments

# prephasing
/restricted/projectnb/cte/_xdhan/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -B /restricted/projectnb/cte/Challenge_Project_2022/data/genetic/${batch}_batch/_mapt_region/bed_files/MAPT.topmed \
        -M /restricted/projectnb/cte/Challenge_Project_2022/data/genetic/reference/mapt/genetic_map_chr17_combined_b37.txt \
        --input-ref /restricted/projectnb/cte/Challenge_Project_2022/data/genetic/reference/mapt/_1000GP_Phase3/1000GP_Phase3_chr17.hap.gz /restricted/projectnb/cte/Challenge_Project_2022/data/genetic/reference/mapt/_1000GP_Phase3/1000GP_Phase3_chr17.legend /restricted/projectnb/cte/Challenge_Project_2022/data/genetic/reference/mapt/_1000GP_Phase3/1000GP_Phase3.sample \
        --exclude-snp /restricted/projectnb/cte/Challenge_Project_2022/imputation/${batch}_batch/prephasing/gwas.alignments.snp.strand.exclude \
        -O /restricted/projectnb/cte/Challenge_Project_2022/imputation/${batch}_batch/prephasing/${batch}_batch.topmed.phased.with.ref \
        -T 28

mv ./shapeit_* /restricted/projectnb/cte/Challenge_Project_2022/imputation/${batch}_batch/prephasing/

echo "==============================================================================================================================================="
echo "Imputation using impute.v2"
echo "==============================================================================================================================================="

mkdir /restricted/projectnb/cte/Challenge_Project_2022/imputation/${batch}_batch/impute

/restricted/projectnb/cte/_xdhan/impute_v2.3.2_x86_64_dynamic/impute2 -use_prephased_g \
        -m /restricted/projectnb/cte/Challenge_Project_2022/data/genetic/reference/mapt/genetic_map_chr17_combined_b37.txt \
        -h /restricted/projectnb/cte/Challenge_Project_2022/data/genetic/reference/mapt/_impute/ref.paper.impute.haps \
        -l /restricted/projectnb/cte/Challenge_Project_2022/data/genetic/reference/mapt/_impute/ref.paper.impute.legend \
        -known_haps_g /restricted/projectnb/cte/Challenge_Project_2022/imputation/${batch}_batch/prephasing/${batch}_batch.topmed.phased.with.ref.haps \
        -int 43865082 45084999 \
        -Ne 20000 \
        -o /restricted/projectnb/cte/Challenge_Project_2022/imputation/${batch}_batch/impute/${batch}_batch.topmed.phased.impute2


echo "==============================================================================================================================================="
echo "Converting imputed data to haplotypic data (2nd phasing)"
echo "==============================================================================================================================================="

mkdir /restricted/projectnb/cte/Challenge_Project_2022/imputation/${batch}_batch/2nd_phasing

echo "convert imputed result into gen format"
Rscript --vanilla /restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20221027/convert.gen.R ${batch}

echo "convert gen format into plink format"
cp /restricted/projectnb/cte/Challenge_Project_2022/imputation/${batch}_batch/prephasing/${batch}_batch.topmed.phased.with.ref.sample /restricted/projectnb/cte/Challenge_Project_2022/imputation/${batch}_batch/2nd_phasing

sed -i 's/0 0 0 -9/NA NA NA NA/g' /restricted/projectnb/cte/Challenge_Project_2022/imputation/${batch}_batch/2nd_phasing/${batch}_batch.topmed.phased.with.ref.sample

plink2 --gen /restricted/projectnb/cte/Challenge_Project_2022/imputation/${batch}_batch/2nd_phasing/${batch}_batch.topmed.phased.impute2 ref-first --make-bed --geno 0.05 --mind 0.05 --out /restricted/projectnb/cte/Challenge_Project_2022/imputation/${batch}_batch/2nd_phasing/${batch}_batch.topmed.phased.impute2 --sample /restricted/projectnb/cte/Challenge_Project_2022/imputation/${batch}_batch/2nd_phasing/${batch}_batch.topmed.phased.with.ref.sample

echo "2nd phasing"
# check
/restricted/projectnb/cte/_xdhan/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -check \
        -B /restricted/projectnb/cte/Challenge_Project_2022/imputation/${batch}_batch/2nd_phasing/${batch}_batch.topmed.phased.impute2 \
        -M /restricted/projectnb/cte/Challenge_Project_2022/data/genetic/reference/mapt/genetic_map_chr17_combined_b37.txt \
        --input-ref /restricted/projectnb/cte/Challenge_Project_2022/data/genetic/reference/mapt/_phasing/imputation_cnv_panel_1kg.reference.convert.hap /restricted/projectnb/cte/Challenge_Project_2022/data/genetic/reference/mapt/_phasing/imputation_cnv_panel_1kg.reference.legend /restricted/projectnb/cte/Challenge_Project_2022/data/genetic/reference/mapt/_phasing/imputation_cnv_panel_1kg.reference.sample \
        --output-log /restricted/projectnb/cte/Challenge_Project_2022/imputation/${batch}_batch/2nd_phasing/gwas.alignments \
# shapeit
/restricted/projectnb/cte/_xdhan/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -B /restricted/projectnb/cte/Challenge_Project_2022/imputation/${batch}_batch/2nd_phasing/${batch}_batch.topmed.phased.impute2 \
        -M /restricted/projectnb/cte/Challenge_Project_2022/data/genetic/reference/mapt/genetic_map_chr17_combined_b37.txt \
        --input-ref /restricted/projectnb/cte/Challenge_Project_2022/data/genetic/reference/mapt/_phasing/imputation_cnv_panel_1kg.reference.convert.hap /restricted/projectnb/cte/Challenge_Project_2022/data/genetic/reference/mapt/_phasing/imputation_cnv_panel_1kg.reference.legend /restricted/projectnb/cte/Challenge_Project_2022/data/genetic/reference/mapt/_phasing/imputation_cnv_panel_1kg.reference.sample \
        --exclude-snp /restricted/projectnb/cte/Challenge_Project_2022/imputation/${batch}_batch/2nd_phasing/gwas.alignments.snp.strand.exclude \
        -O /restricted/projectnb/cte/Challenge_Project_2022/imputation/${batch}_batch/2nd_phasing/${batch}_batch.topmed.phased.impute2.phased.with.ref \
        -T 28

mv shapeit_* /restricted/projectnb/cte/Challenge_Project_2022/imputation/${batch}_batch/2nd_phasing/


echo "==============================================================================================================================================="
echo "Generating the subhaplotypes"
echo "==============================================================================================================================================="

mkdir /restricted/projectnb/cte/Challenge_Project_2022/imputation/${batch}_batch/subhaplotypes

grep "chr17:" /restricted/projectnb/cte/Challenge_Project_2022/imputation/${batch}_batch/2nd_phasing/${batch}_batch.topmed.phased.impute2.phased.with.ref.haps > /restricted/projectnb/cte/Challenge_Project_2022/imputation/${batch}_batch/subhaplotypes/${batch}_batch.topmed.phased.impute2.phased.with.ref.haps.surrogate

grep "chr17:" /restricted/projectnb/cte/Challenge_Project_2022/imputation/${batch}_batch/impute/${batch}_batch.topmed.phased.impute2 > /restricted/projectnb/cte/Challenge_Project_2022/imputation/${batch}_batch/subhaplotypes/${batch}_batch.topmed.phased.impute2.surrogate.dosage
echo "convert haplotypic to subhaplotypic data"
Rscript --vanilla /restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20221027/subhaplotypes.R ${batch}

echo "Finished All!"
