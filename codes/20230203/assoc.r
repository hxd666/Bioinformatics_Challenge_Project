# association test using subhaplotypes
rm(list = ls())
require(mice)
require(MASS)
require(rms)

# ==============================================================================================
subject <- 'all'

if (!subject %in% c('all', 'cases', 'age')){
  stop('subject (1st) parameter should be either age, cases, or age!', )
}
# ==============================================================================================

# Load Data
load('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20230201/imputatoin.RData')
# Set PCs formula
pred.list   <- read.table('predictors_list', as.is = T)
pred.list   <- pred.list$V1
pcs.formula <- sapply(1:10, function(x){paste0('PC',x)})
pcs.formula <- paste(pcs.formula, collapse = ' + ')
aa_subset <- pheno$race %in% 2

# ==============================================================================================
model.all.data <- data.frame('Subjects' = 0, 'Model' = 0, 'Phenotype' = 0, 'Predictor'= 0, 'Freq' = 0,
                             'N' = 0, 'Effect' = 0, 'SE' = 0, 'P-value' = 0, stringsAsFactors = F)
# All Subjects Models
for (predictor in c('H1B1G1')){
  
  for (.pheno in c(pheno_list_binary, pheno_list_categorical, pheno_list_continuous)){
    
    print(paste('Model All:', .pheno))
    # Model 1: .phenotype ~ Predictor + Age + PCs
    model1.form <- paste(.pheno,  " ~ ",  predictor, " + agedeath + ", pcs.formula, sep = '')
    
    # Fit One:
    subset <- !is.na(pheno[, predictor]) & aa_subset
    if (.pheno %in% c("CTE", "DementiaHx")){
      fit <- with(imp, glm(as.formula(model1.form), na.action = "na.omit", family = binomial(link = "logit"), subset = subset))
    } else if (.pheno %in% pheno_list_categorical){
      fit <- with(imp, polr(as.formula(model1.form), na.action = "na.omit", Hess = T, subset = subset))
    } else if (.pheno %in% pheno_list_continuous){
      fit <- with(imp, glm(as.formula(model1.form), na.action = "na.omit", family = gaussian(link = "identity"), subset = subset))
    }
    pool.fit <- pool(fit)
    result <- summary(pool.fit)
    assoc_subject <- !(is.na(complete(imp,1)[,.pheno])) & subset
    n <- sum(assoc_subject)
    effect <- round(result$estimate[result$term == predictor], 3)
    se <- round(result$std.error[result$term == predictor], 3)
    p  <- result$p.value[result$term == predictor]
    freq <- round(sum(pheno[,predictor][assoc_subject])/(2*n), 3)
    output.row <- c('African', '1', .pheno, predictor, freq, n, effect, se, p)
    output.row <- unname(output.row)
    model.all.data <- rbind(model.all.data, output.row)
  }
}

# Calc FDR
model.all.data <- model.all.data[2:nrow(model.all.data), ]
model.all.data$Freq >= 0.05
model.all.data$FDR <- p.adjust(model.all.data$P.value, method = 'fdr')

# save result
outputname <- 'results/Model1.AFR.assoc.txt'
write.table(model.all.data, outputname, col.names = T, row.names = F, sep = '\t', quote = F)
# ==============================================================================================


