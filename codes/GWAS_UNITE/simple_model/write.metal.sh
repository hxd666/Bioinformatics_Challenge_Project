#!/bin/bash -l
risk=${1} # risk should be CTE or CTEstage

if [ -d metal ]; then
	echo 'metal/ already existed!'
else
	mkdir metal
fi

if [ -d metal/${risk} ]; then
	echo "metal/${risk}/ already existed!" 
else
	mkdir metal/${risk}
fi

echo '#!/bin/bash -l' > metal/${risk}/metal.input

echo "# run metal
GENOMICCONTROL ON
AVERAGEFREQ ON
scheme STDERR
CUSTOMVARIABLE N
#CUSTOMVARIABLE N.case
#CUSTOMVARIABLE N.ctrl

OUTFILE UNITE_${risk}.assoc. .txt

ADDFILTER BETA.snp >= -5
ADDFILTER BETA.snp <= 5
ADDFILTER FREQ >= 0.05
ADDFILTER FREQ <= 0.95
MARKER LocID
ALLELE A1 A2
FREQ FREQ
PVALUE P.snp
EFFECT BETA.snp
STDERR SE.snp

PROCESS /restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/GWAS_UNITE/endophenotypes/batch1/${risk}/UNITE_batch1_${risk}.assoc.gz
PROCESS /restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/GWAS_UNITE/endophenotypes/batch2/${risk}/UNITE_batch2_${risk}.assoc.gz

ANALYZE HETEROGENEITY" >> metal/${risk}/metal.input

echo '#!/bin/bash -l
module load metal
metal < metal.input' > metal/${risk}/run.metal.sh

chmod a+x metal/${risk}/run.metal.sh

cp /restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/GWAS_UNITE/simple_model/metal/CTE/submit.run.filter_metalout_df.v0.sh metal/${risk}/
