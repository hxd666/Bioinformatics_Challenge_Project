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
load('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20230515/imputation.RData')
# Set PCs formula
pred.list   <- 'rs6910507'
pcs.formula <- sapply(1:10, function(x){paste0('PC',x)})
pcs.formula <- paste(pcs.formula, collapse = ' + ')
euro_subset <- pheno$race %in% 1 & pheno$Batch == "I"

# ==============================================================================================
model.all.data <- data.frame('Subjects' = 0, 'Model' = 0, 'Phenotype' = 0, 'Predictor'= 0, 'Freq' = 0,
                             'N' = 0, 'Effect' = 0, 'SE' = 0, 'P-value' = 0, stringsAsFactors = F)
# All Subjects Models
for (predictor in pred.list){
  
  for (.pheno in c(pheno_list_binary, pheno_list_categorical, pheno_list_continuous)){
    
    print(paste('Model All:', .pheno))
    # Model 1: .phenotype ~ Predictor + Age + PCs
    model1.form <- paste(.pheno,  " ~ ",  predictor, " + agedeath + ", pcs.formula, sep = '')
    
    # Fit One:
    subset <- !is.na(pheno[, predictor]) & euro_subset
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
    effect <- round(result$estimate[rownames(result) == predictor], 3)
    se <- round(result$std.error[rownames(result) == predictor], 3)
    p  <- result$p.value[rownames(result) == predictor]
    freq <- round(sum(pheno[,predictor][assoc_subject])/(2*n), 3)
    output.row <- c('Euro', '1', .pheno, predictor, freq, n, effect, se, p)
    output.row <- unname(output.row)
    model.all.data <- rbind(model.all.data, output.row)
  }
}

# save
model.all.data <- model.all.data[2:nrow(model.all.data), ]
model.all.data$Batch = 'I'
outputname <- 'results/Model1.Euro.batchI.assoc.txt'
write.table(model.all.data, outputname, col.names = T, row.names = F, sep = '\t', quote = F)






euro_subset <- pheno$race %in% 1 & pheno$Batch == "II"
# ==============================================================================================
model.all.data <- data.frame('Subjects' = 0, 'Model' = 0, 'Phenotype' = 0, 'Predictor'= 0, 'Freq' = 0,
                             'N' = 0, 'Effect' = 0, 'SE' = 0, 'P-value' = 0, stringsAsFactors = F)
# All Subjects Models
for (predictor in pred.list){
  
  for (.pheno in c(pheno_list_binary, pheno_list_categorical, pheno_list_continuous)){
    
    print(paste('Model All:', .pheno))
    # Model 1: .phenotype ~ Predictor + Age + PCs
    model1.form <- paste(.pheno,  " ~ ",  predictor, " + agedeath + ", pcs.formula, sep = '')
    
    # Fit One:
    subset <- !is.na(pheno[, predictor]) & euro_subset
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
    effect <- round(result$estimate[rownames(result) == predictor], 3)
    se <- round(result$std.error[rownames(result) == predictor], 3)
    p  <- result$p.value[rownames(result) == predictor]
    freq <- round(sum(pheno[,predictor][assoc_subject])/(2*n), 3)
    output.row <- c('Euro', '1', .pheno, predictor, freq, n, effect, se, p)
    output.row <- unname(output.row)
    model.all.data <- rbind(model.all.data, output.row)
  }
}
model.all.data <- model.all.data[2:nrow(model.all.data), ]
model.all.data$Batch = 'II'
outputname <- 'results/Model1.Euro.batchII.assoc.txt'
write.table(model.all.data, outputname, col.names = T, row.names = F, sep = '\t', quote = F)
# ==============================================================================================




euro_subset <- pheno$race %in% 1 & pheno$Batch == "I" & pheno$CTE %in% c(1)

# ==============================================================================================
model.all.data <- data.frame('Subjects' = 0, 'Model' = 0, 'Phenotype' = 0, 'Predictor'= 0, 'Freq' = 0,
                             'N' = 0, 'Effect' = 0, 'SE' = 0, 'P-value' = 0, stringsAsFactors = F)
# All Subjects Models
for (predictor in pred.list){
  
  for (.pheno in c(pheno_list_binary, pheno_list_categorical, pheno_list_continuous)){
    
    print(paste('Model All:', .pheno))
    # Model 1: .phenotype ~ Predictor + Age + PCs
    model1.form <- paste(.pheno,  " ~ ",  predictor, " + agedeath + ", pcs.formula, sep = '')
    
    # Fit One:
    subset <- !is.na(pheno[, predictor]) & euro_subset
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
    effect <- round(result$estimate[rownames(result) == predictor], 3)
    se <- round(result$std.error[rownames(result) == predictor], 3)
    p  <- result$p.value[rownames(result) == predictor]
    freq <- round(sum(pheno[,predictor][assoc_subject])/(2*n), 3)
    output.row <- c('Euro', '1', .pheno, predictor, freq, n, effect, se, p)
    output.row <- unname(output.row)
    model.all.data <- rbind(model.all.data, output.row)
  }
}

# save
model.all.data <- model.all.data[2:nrow(model.all.data), ]
model.all.data$Batch = 'I'
outputname <- 'results/Model1.Euro.batchI.CTEonly.assoc.txt'
write.table(model.all.data, outputname, col.names = T, row.names = F, sep = '\t', quote = F)






euro_subset <- pheno$race %in% 1 & pheno$Batch == "II" & pheno$CTE %in% c(1)
# ==============================================================================================
model.all.data <- data.frame('Subjects' = 0, 'Model' = 0, 'Phenotype' = 0, 'Predictor'= 0, 'Freq' = 0,
                             'N' = 0, 'Effect' = 0, 'SE' = 0, 'P-value' = 0, stringsAsFactors = F)
# All Subjects Models
for (predictor in pred.list){
  
  for (.pheno in c(pheno_list_binary, pheno_list_categorical, pheno_list_continuous)){
    
    print(paste('Model All:', .pheno))
    # Model 1: .phenotype ~ Predictor + Age + PCs
    model1.form <- paste(.pheno,  " ~ ",  predictor, " + agedeath + ", pcs.formula, sep = '')
    
    # Fit One:
    subset <- !is.na(pheno[, predictor]) & euro_subset
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
    effect <- round(result$estimate[rownames(result) == predictor], 3)
    se <- round(result$std.error[rownames(result) == predictor], 3)
    p  <- result$p.value[rownames(result) == predictor]
    freq <- round(sum(pheno[,predictor][assoc_subject])/(2*n), 3)
    output.row <- c('Euro', '1', .pheno, predictor, freq, n, effect, se, p)
    output.row <- unname(output.row)
    model.all.data <- rbind(model.all.data, output.row)
  }
}
model.all.data <- model.all.data[2:nrow(model.all.data), ]
model.all.data$Batch = 'II'
outputname <- 'results/Model1.Euro.batchII.CTEonly.assoc.txt'
write.table(model.all.data, outputname, col.names = T, row.names = F, sep = '\t', quote = F)
