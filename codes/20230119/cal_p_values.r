# calculate the p-value for permutation tests in 20230116
setwd('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20230119')

predictor <- commandArgs(T)[1]
message(predictor)

# Load Data
stats_ds <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20221211_2/results/Model1.Euro.assoc.txt',
                       header = T, as.is = T)
stats_ds <- subset(stats_ds, Predictor == predictor)

pheno <- stats_ds$Phenotype

p_values_fn <- '/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20230116/results/XXXX.permutation.P.values'
p_values_fn <- gsub('XXXX', predictor, p_values_fn)
p_values_ds <- read.table(p_values_fn, header = T, as.is = T)

p_values_ds$min.P <- apply(p_values_ds[, pheno], 1, min)
  
# Calc P-values for Permutation tests
stats_ds$Perm_P <- NA

for (.pheno in pheno){
  .per.p  <- mean(p_values_ds$min.P < stats_ds$P.value[stats_ds$Phenotype == .pheno], na.rm = T)
  stats_ds$Perm_P[stats_ds$Phenotype == .pheno] <- .per.p
}

output_name <- '/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20230119/results/XXXX.Model1.Euro.permutation.assoc.txt'
output_name <- gsub('XXXX', predictor, output_name)
write.table(stats_ds, output_name, col.names = T, row.names = F, quote = F, sep = '\t')
