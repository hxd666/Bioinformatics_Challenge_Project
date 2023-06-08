#!/bin/bash -l

module load R/4.0.0

Rscript --vanilla permutation.batchII.r $1
