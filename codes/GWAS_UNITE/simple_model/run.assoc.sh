#!/bin/bash -l
risk=$1 # risk should be CTE or CTEstage
study=$2 # study should be batch1 or batch2
chr=$3 # chr should be any integer from 1 to 22

module load htslib/1.16 bcftools/1.16 R/4.0.0

Rscript --vanilla assoc.R $risk $study $chr
