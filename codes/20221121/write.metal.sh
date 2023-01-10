#!/bin/bash -l

echo "# run metal
GENOMICCONTROL OFF
AVERAGEFREQ OFF
scheme STDERR
CUSTOMVARIABLE N

OUTFILE output/model${1}.euro.metal. .txt

#ADDFILTER BETA.int >= -5
#ADDFILTER BETA.int <= 5
#ADDFILTER FREQ >= 0.05
#ADDFILTER FREQ <= 0.95
MARKER Phenotype
#ALLELE Phenotype Predictor
#FREQ FREQ
PVALUE P.value
EFFECT Effect
STDERR SE
PROCESS input/model${1}.batchI.euro.txt
PROCESS input/model${1}.batchII.euro.txt
ANALYZE HETEROGENEITY" > /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_assoc/20221121/metal.model${1}.input

head -n 1 /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_assoc/20221118/model${1}.euro.txt > /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_assoc/20221121/input/model${1}.batchI.euro.txt
grep -w 'I' /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_assoc/20221118/model${1}.euro.txt >> /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_assoc/20221121/input/model${1}.batchI.euro.txt

head -n 1 /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_assoc/20221118/model${1}.euro.txt > /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_assoc/20221121/input/model${1}.batchII.euro.txt
grep -w 'II' /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_assoc/20221118/model${1}.euro.txt >> /restricted/projectnb/cte/Challenge_Project_2022/xdhan/_assoc/20221121/input/model${1}.batchII.euro.txt
