#!/bin/bash -l

module load R/4.0.0

Rscript --vanilla permutation.batchI.r $1
