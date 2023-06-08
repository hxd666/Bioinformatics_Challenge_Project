# association test using subhaplotypes
rm(list = ls())
require(mice)
require(MASS)
require(rms)

predictor <- commandArgs(T)[1]
# ==============================================================================================
subject <- 'all'
set.seed(100)

if (!subject %in% c('all', 'cases', 'age')){
  stop('subject (1st) parameter should be either age, cases, or age!', )
}
# ==============================================================================================

# Load Data
load('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20230201/imputatoin.RData')
# Set PCs formula
pcs.formula <- sapply(1:10, function(x){paste0('PC',x)})
pcs.formula <- paste(pcs.formula, collapse = ' + ')
euro_subset <- pheno$race %in% 1

# ==============================================================================================
# Permutation Test
L <- 1000
N <- nrow(complete(imp, 1)[complete(imp,1)$race %in% "1" & complete(imp,1)$Batch %in% 'I', ])
model.all.simulation <- data.frame('Subjects' = 0, 'Model' = 0, 'SNP' = 0, 'Freq' = 0,
                                   'N' = 0, 'Effect' =0, 'SE' = 0, 'P.value' = 0, stringsAsFactors = F)

for (l in (1:L)) {
  cat(l); cat('.')
  l.sample <- sample(N)
  
  assoc.list <- list()
  for (iter in (1:10)){
    l.imp.iter <- complete(imp, iter)[complete(imp,1)$race %in% "1" & complete(imp,1)$Batch %in% 'I', ]
    l.imp.iter.outcomes <- l.imp.iter[l.sample, c(pheno_list_binary, pheno_list_categorical, pheno_list_continuous)]
    l.imp.iter.others   <- l.imp.iter[, !colnames(l.imp.iter) %in% c(pheno_list_binary, pheno_list_categorical, pheno_list_continuous)]
    l.imp.iter.assoc    <- cbind(l.imp.iter.outcomes, l.imp.iter.others)
    assoc.list <- append(assoc.list, list(l.imp.iter.assoc))
  }
  
  #p.vec <- c()
  for (.pheno in c(pheno_list_binary, pheno_list_categorical, pheno_list_continuous)){
    model1.form <- paste(.pheno,  " ~ ",  predictor, " + agedeath + ", pcs.formula, sep = '')
    subset <- !is.na(pheno[, predictor]) & euro_subset & pheno$Batch %in% 'I'
    if (.pheno %in% c("CTE", "DementiaHx")){
      fit <- lapply(assoc.list, function(x){glm(as.formula(model1.form), x, na.action="na.omit", family=binomial(link="logit"))})
    } else if (.pheno %in% pheno_list_categorical){
      fit <- lapply(assoc.list, function(x){polr(as.formula(model1.form), x, na.action="na.omit", Hess = T)})
    } else if (.pheno %in% pheno_list_continuous){
      fit <- lapply(assoc.list, function(x){glm(as.formula(model1.form), x, na.action="na.omit", family=gaussian(link="identity"))})
      }
    pool.fit <- pool(fit)
    result <- summary(pool.fit)
    snp <- paste(.pheno, predictor, l, sep = '_')
    n   <- sum(subset)
    freq <- sum(pheno[subset, ]$H1B1G1)/(2*n)
    effect <- result$estimate[result$term == predictor]
    se     <- result$std.error[result$term == predictor]
    p      <- result$p.value[result$term == predictor]
    output.row <- c('Euro', '1', snp, freq, n, effect, se, p)
    model.all.simulation <- rbind(model.all.simulation, output.row)
    #p.vec <- c(p.vec, p)
  }
  #output.row <- c(predictor, l, p.vec)  
  #model.all.simulation <- rbind(model.all.simulation, output.row)
}
model.all.simulation <- model.all.simulation[2:nrow(model.all.simulation), ]
model.all.simulation$A1 <- 'A'
model.all.simulation$A2 <- 'G'
# ==============================================================================================

output.filename <- 'permutations/XXXX.permutation.batch1.P.values'
output.filename <- gsub('XXXX', predictor, output.filename)
write.table(model.all.simulation, output.filename, row.names = F, col.names = T, quote = F, 
            sep = '\t')
