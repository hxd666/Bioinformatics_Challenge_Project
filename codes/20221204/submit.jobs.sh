#!/bin/bash -l

module load R/4.0.0

while read -r predictor; do
	echo "Predictor: $predictor"
	Rscript --vanilla assoc.r all $predictor
	Rscript --vanilla assoc.r cases $predictor
	Rscript --vanilla assoc.r age $predictor
done < predictors_list
