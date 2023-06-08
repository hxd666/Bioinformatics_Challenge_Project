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
load('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20221128/imputatoin.RData')
# Set PCs formula
pcs.formula <- sapply(1:10, function(x){paste0('PC',x)})
pcs.formula <- paste(pcs.formula, collapse = ' + ')
euro_subset <- pheno$race %in% 1

# ==============================================================================================
# Permutation Test
L <- 1000
n <- nrow(complete(imp, 1))
model.all.simulation <- matrix(0, nrow = 1, ncol = (2+length(c(pheno_list_binary, pheno_list_categorical))))
model.all.simulation <- as.data.frame(model.all.simulation)
colnames(model.all.simulation) <- c('Predictor', 'Simulation', pheno_list_binary, pheno_list_categorical)

for (l in (1:L)) {
  cat(l); cat('.')
  l.sample <- sample(n)
  
  assoc.list <- list()
  for (iter in (1:10)){
    l.imp.iter <- complete(imp, iter)
    l.imp.iter.outcomes <- l.imp.iter[l.sample, c(pheno_list_binary, pheno_list_categorical)]
    l.imp.iter.others   <- l.imp.iter[, !colnames(l.imp.iter) %in% c(pheno_list_binary, pheno_list_categorical)]
    l.imp.iter.assoc    <- cbind(l.imp.iter.outcomes, l.imp.iter.others)
    assoc.list <- append(assoc.list, list(l.imp.iter.assoc))
  }
  
  p.vec <- c()
  for (.pheno in c(pheno_list_binary, pheno_list_categorical)){
    model1.form <- paste(.pheno,  " ~ ",  predictor, " + agedeath + ", pcs.formula, sep = '')
    subset <- !is.na(pheno[, predictor]) & euro_subset
    if (.pheno %in% c("CTE", "DementiaHx")){
      fit <- lapply(assoc.list, function(x){glm(as.formula(model1.form), x, na.action="na.omit", family=binomial(link="logit"), subset = subset)})
    } else {
      fit <- lapply(assoc.list, function(x){polr(as.formula(model1.form), x, na.action="na.omit", Hess = T, subset = subset)})
    }
    pool.fit <- pool(fit)
    result <- summary(pool.fit)
    p  <- result$p.value[result$term == predictor]
    p.vec <- c(p.vec, p)
  }
  output.row <- c(predictor, l, p.vec)  
  model.all.simulation <- rbind(model.all.simulation, output.row)
}
model.all.simulation <- model.all.simulation[2:nrow(model.all.simulation), ]
# ==============================================================================================

output.filename <- 'results/XXXX.permutation.P.values'
output.filename <- gsub('XXXX', predictor, output.filename)
write.table(model.all.simulation, output.filename, row.names = F, col.names = T, quote = F, 
            sep = '\t')
