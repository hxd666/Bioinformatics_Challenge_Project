# calculate the p-value for permutation tests in 20230116
setwd('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20230201')

predictor <- commandArgs(T)[1]
message(predictor)

# Load Data
stats_ds <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20230201/results/Model1.Euro.assoc.txt',
                       header = T, as.is = T)
stats_ds <- subset(stats_ds, Predictor == predictor)

pheno <- stats_ds$Phenotype
pheno_global <- pheno[c(1, 2, 13, 15)]
pheno_regional <- pheno[! pheno %in% pheno_global]

p_values_fn <- '/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20230201/permutations/XXXX.permutation.P.values'
p_values_fn <- gsub('XXXX', predictor, p_values_fn)
p_values_ds <- read.table(p_values_fn, header = T, as.is = T)

p_values_ds$min.P.global <- apply(p_values_ds[, pheno_global], 1, min)
p_values_ds$min.P.regional <- apply(p_values_ds[, pheno_regional], 1, min)

# Calc P-values for Permutation tests
stats_ds$Perm_P <- NA

for (.pheno in pheno){
  if (.pheno %in% pheno_global){
    .per.p  <- mean(p_values_ds$min.P.global < stats_ds$P.value[stats_ds$Phenotype == .pheno], na.rm = T)
  } else if (.pheno %in% pheno_regional){
    .per.p  <- mean(p_values_ds$min.P.regional < stats_ds$P.value[stats_ds$Phenotype == .pheno], na.rm = T)
  }
  stats_ds$Perm_P[stats_ds$Phenotype == .pheno] <- .per.p
}

output_name <- '/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20230201/p_adjust/XXXX.Model1.Euro.permutation.assoc.txt'
output_name <- gsub('XXXX', predictor, output_name)
write.table(stats_ds, output_name, col.names = T, row.names = F, quote = F, sep = '\t')
