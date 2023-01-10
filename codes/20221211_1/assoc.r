# association test using subhaplotypes
require(mice)
require(MASS)
require(rms)

# ==============================================================================================
# Set Arguments
t.args <- commandArgs(T)
subject <- t.args[1] 
# subject <- 'cases'
predictor <- t.args[2]
#predictor <- 'H1B1G1'

if (!subject %in% c('all', 'cases', 'age')){
  stop('subject (1st) parameter should be either age, cases, or age!', )
}
# ==============================================================================================

# Load Data
load('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20221128/imputatoin.RData')
# Set PCs formula
pcs.formula <- sapply(1:10, function(x){paste0('PC',x)})
pcs.formula <- paste(pcs.formula, collapse = ' + ')

european_subset <- pheno$race %in% 1

# ==============================================================================================
# All Subjects Models
if (subject == 'all'){
  model.all.data <- data.frame('Subjects' = 0, 'Model' = 0, 'Phenotype' = 0, 'Predictor'= 0, 'Freq' = 0,
                               'N' = 0, 'Effect' = 0, 'SE' = 0, 'P-value' = 0, stringsAsFactors = F)
  
  for (.pheno in pheno_list_categorical){
    print(paste('Model All:', .pheno))
    # Model 1: .phenotype ~ Predictor + Age + PCs
    model1.form <- paste(.pheno,  " ~ ",  predictor, " + agedeath + ", pcs.formula, sep = '')
    # Model 2: .phenotype ~ Predictor + Age + PCs + Duration (Football Only)
    model2.form <- paste(.pheno,  " ~ ",  predictor, " + agedeath + ", pcs.formula, ' + footyrs', sep = '')
    # Model 3: .phenotype ~ Predictor + Age + PCs + Duration + Duration*Predictor (Football Only)
    model3.form <- paste(.pheno,  " ~ ",  predictor, " + agedeath + ", pcs.formula, ' + footyrs + footyrs*', predictor,sep = '')
    
    # Fit One:
    subset <- !is.na(pheno[, predictor]) & european_subset
    fit <- with(imp, polr(as.formula(model1.form), na.action="na.omit", Hess = T, subset = subset))
    pool.fit <- pool(fit)
    result <- summary(pool.fit)
    assoc_subject <- !(is.na(complete(imp,1)[,.pheno])) & subset
    n <- sum(assoc_subject)
    effect <- round(result$estimate[result$term == predictor], 3)
    se <- round(result$std.error[result$term == predictor], 3)
    p  <- result$p.value[result$term == predictor]
    freq <- round(sum(pheno[,predictor][assoc_subject])/(2*n), 3)
    output.row <- c('All', '1', .pheno, predictor, freq, n, effect, se, p)
    output.row <- unname(output.row)
    model.all.data <- rbind(model.all.data, output.row)
    
    # Fit Two:
    football_subset <- !is.na(pheno$footyrs) & european_subset
    subset <- subset & football_subset
    fit <- with(imp, polr(as.formula(model2.form), na.action="na.omit", Hess = T, subset = subset))
    pool.fit <- pool(fit)
    result <- summary(pool.fit)
    assoc_subject <- !(is.na(complete(imp,1)[,.pheno])) & subset
    n <- sum(assoc_subject)
    effect <- round(result$estimate[result$term == predictor], 3)
    se <- round(result$std.error[result$term == predictor], 3)
    p  <- result$p.value[result$term == predictor]
    freq <- round(sum(pheno[,predictor][assoc_subject])/(2*n), 3)
    output.row <- c('All', '2', .pheno, predictor, freq, n, effect, se, p)
    output.row <- unname(output.row)
    model.all.data <- rbind(model.all.data, output.row)
    
    # Fit Three:
    fit <- with(imp, polr(as.formula(model3.form), na.action="na.omit", Hess = T, subset = subset))
    pool.fit <- pool(fit)
    result <- summary(pool.fit)
    assoc_subject <- !(is.na(complete(imp,1)[,.pheno])) & subset
    n <- sum(assoc_subject)
    effect <- round(result$estimate[result$term == paste0(predictor, ":footyrs")], 3)
    se <- round(result$std.error[result$term == paste0(predictor, ":footyrs")], 3)
    p  <- result$p.value[result$term == paste0(predictor, ":footyrs")]
    freq <- round(sum(pheno[,predictor][assoc_subject])/(2*n), 3)
    output.row <- c('All', '3', .pheno, paste0(predictor, "*footyrs"), freq, n, effect, se, p)
    output.row <- unname(output.row)
    model.all.data <- rbind(model.all.data, output.row)
  }
  # Calc FDR
  model.all.data <- model.all.data[2:nrow(model.all.data), ]
  model.all.data$FDR <- NA
  model.all.data$FDR[model.all.data$Model == '1'] <- p.adjust(model.all.data$P.value[model.all.data$Model == '1'], method = 'fdr')
  model.all.data$FDR[model.all.data$Model == '2'] <- p.adjust(model.all.data$P.value[model.all.data$Model == '2'], method = 'fdr')
  model.all.data$FDR[model.all.data$Model == '3'] <- p.adjust(model.all.data$P.value[model.all.data$Model == '3'], method = 'fdr')
  # save result
  outputname <- 'results/Model.All.XXXX.assoc.txt'
  outputname <- gsub('XXXX', predictor, outputname)
  write.table(model.all.data, outputname, col.names = T, row.names = F, sep = '\t', quote = F)
}
# ==============================================================================================


# ==============================================================================================
# Cases Subjects Models
if (subject == 'cases'){
  model.cases.data <- data.frame('Subjects' = 0, 'Model' = 0, 'Phenotype' = 0, 'Predictor'= 0, 'Freq' = 0,
                               'N' = 0, 'Effect' = 0, 'SE' = 0, 'P-value' = 0, stringsAsFactors = F)
  
  for (.pheno in pheno_list_categorical){
    print(paste('Model Cases:', .pheno))
    # Model 1: .phenotype ~ Predictor + Age + PCs
    model1.form <- paste(.pheno,  " ~ ",  predictor, " + agedeath + ", pcs.formula, sep = '')
    # Model 2: .phenotype ~ Predictor + Age + PCs + Duration (Footbcases Only)
    model2.form <- paste(.pheno,  " ~ ",  predictor, " + agedeath + ", pcs.formula, ' + footyrs', sep = '')
    # Model 3: .phenotype ~ Predictor + Age + PCs + Duration + Duration*Predictor (Footbcases Only)
    model3.form <- paste(.pheno,  " ~ ",  predictor, " + agedeath + ", pcs.formula, ' + footyrs + footyrs*', predictor,sep = '')
    
    # Fit One:
    subset <- !is.na(pheno[, predictor]) & pheno$CTE == '1'  & european_subset
    fit <- with(imp, polr(as.formula(model1.form), na.action="na.omit", Hess = T, subset = subset))
    pool.fit <- pool(fit)
    result <- summary(pool.fit)
    assoc_subject <- !(is.na(complete(imp,1)[,.pheno])) & subset
    n <- sum(assoc_subject)
    effect <- round(result$estimate[result$term == predictor], 3)
    se <- round(result$std.error[result$term == predictor], 3)
    p  <- result$p.value[result$term == predictor]
    freq <- round(sum(pheno[,predictor][assoc_subject])/(2*n), 3)
    output.row <- c('Cases', '1', .pheno, predictor, freq, n, effect, se, p)
    output.row <- unname(output.row)
    model.cases.data <- rbind(model.cases.data, output.row)
    
    # Fit Two:
    footbcases_subset <- !is.na(pheno$footyrs)  & european_subset
    subset <- subset & footbcases_subset
    fit <- with(imp, polr(as.formula(model2.form), na.action="na.omit", Hess = T, subset = subset))
    pool.fit <- pool(fit)
    result <- summary(pool.fit)
    assoc_subject <- !(is.na(complete(imp,1)[,.pheno])) & subset
    n <- sum(assoc_subject)
    effect <- round(result$estimate[result$term == predictor], 3)
    se <- round(result$std.error[result$term == predictor], 3)
    p  <- result$p.value[result$term == predictor]
    freq <- round(sum(pheno[,predictor][assoc_subject])/(2*n), 3)
    output.row <- c('Cases', '2', .pheno, predictor, freq, n, effect, se, p)
    output.row <- unname(output.row)
    model.cases.data <- rbind(model.cases.data, output.row)
    
    # Fit Three:
    fit <- with(imp, polr(as.formula(model3.form), na.action="na.omit", Hess = T, subset = subset))
    pool.fit <- pool(fit)
    result <- summary(pool.fit)
    assoc_subject <- !(is.na(complete(imp,1)[,.pheno])) & subset
    n <- sum(assoc_subject)
    effect <- round(result$estimate[result$term == paste0(predictor, ":footyrs")], 3)
    se <- round(result$std.error[result$term == paste0(predictor, ":footyrs")], 3)
    p  <- result$p.value[result$term == paste0(predictor, ":footyrs")]
    freq <- round(sum(pheno[,predictor][assoc_subject])/(2*n), 3)
    output.row <- c('Cases', '3', .pheno, paste0(predictor, "*footyrs"), freq, n, effect, se, p)
    output.row <- unname(output.row)
    model.cases.data <- rbind(model.cases.data, output.row)
  }
  # Calc FDR
  model.cases.data <- model.cases.data[2:nrow(model.cases.data), ]
  model.cases.data$FDR <- NA
  model.cases.data$FDR[model.cases.data$Model == '1'] <- p.adjust(model.cases.data$P.value[model.cases.data$Model == '1'], method = 'fdr')
  model.cases.data$FDR[model.cases.data$Model == '2'] <- p.adjust(model.cases.data$P.value[model.cases.data$Model == '2'], method = 'fdr')
  model.cases.data$FDR[model.cases.data$Model == '3'] <- p.adjust(model.cases.data$P.value[model.cases.data$Model == '3'], method = 'fdr')
  # save result
  outputname <- 'results/Model.Cases.XXXX.assoc.txt'
  outputname <- gsub('XXXX', predictor, outputname)
  write.table(model.cases.data, outputname, col.names = T, row.names = F, sep = '\t', quote = F)
}
# ==============================================================================================


# ==============================================================================================
# Age-Stratified Models
if (subject == 'age'){
  model.age.data <- data.frame('Subjects' = 0, 'Age' = 0, 'Model' = 0, 'Phenotype' = 0, 'Predictor'= 0, 'Freq' = 0,
                               'N' = 0, 'Effect' = 0, 'SE' = 0, 'P-value' = 0, stringsAsFactors = F)
  for (.pheno in pheno_list_categorical){
    print(paste('Model Age-Stratified:', .pheno))
    # Model 1: .phenotype ~ Predictor + Age + PCs
    model1.form <- paste(.pheno,  " ~ ",  predictor, " + agedeath + ", pcs.formula, sep = '')
    # Model 2: .phenotype ~ Predictor + Age + PCs + Duration (Football Only)
    model2.form <- paste(.pheno,  " ~ ",  predictor, " + agedeath + ", pcs.formula, ' + footyrs', sep = '')
    # Model 3: .phenotype ~ Predictor + Age + PCs + Duration + Duration*Predictor (Footbcases Only)
    model3.form <- paste(.pheno,  " ~ ",  predictor, " + agedeath + ", pcs.formula, ' + footyrs + footyrs*', predictor,sep = '')
    # Model 4: .phenotype ~ Predictor + Age_Strata + Predictor*Age_Strata + PCs
    model4.form <- paste(.pheno,  " ~ ", predictor, " + stratum + ", pcs.formula, ' + stratum*', predictor, sep = '')

    # Fit One All:
    subset_older <- !is.na(pheno[, predictor]) & pheno$agedeath > median_age  & european_subset
    subset_younger <- !is.na(pheno[, predictor]) & pheno$agedeath <= median_age  & european_subset
    # ---- older:
    fit <- with(imp, polr(as.formula(model1.form), na.action="na.omit", Hess = T, subset = subset_older))
    pool.fit <- pool(fit)
    result <- summary(pool.fit)
    assoc_subject <- !(is.na(complete(imp,1)[,.pheno])) & subset_older
    n <- sum(assoc_subject)
    effect <- round(result$estimate[result$term == predictor], 3)
    se <- round(result$std.error[result$term == predictor], 3)
    p  <- result$p.value[result$term == predictor]
    freq <- round(sum(pheno[,predictor][assoc_subject])/(2*n), 3)
    output.row <- c('All', 'Older', '1', .pheno, predictor, freq, n, effect, se, p)
    output.row <- unname(output.row)
    model.age.data <- rbind(model.age.data, output.row)
    # ---- younger:
    fit <- with(imp, polr(as.formula(model1.form), na.action="na.omit", Hess = T, subset = subset_younger))
    pool.fit <- pool(fit)
    result <- summary(pool.fit)
    assoc_subject <- !(is.na(complete(imp,1)[,.pheno])) & subset_younger
    n <- sum(assoc_subject)
    effect <- round(result$estimate[result$term == predictor], 3)
    se <- round(result$std.error[result$term == predictor], 3)
    p  <- result$p.value[result$term == predictor]
    freq <- round(sum(pheno[,predictor][assoc_subject])/(2*n), 3)
    output.row <- c('All', 'Younger', '1', .pheno, predictor, freq, n, effect, se, p)
    output.row <- unname(output.row)
    model.age.data <- rbind(model.age.data, output.row)
    # Fit One Cases:
    subset_older <- !is.na(pheno[, predictor]) & pheno$agedeath > median_age & pheno$CTE == '1'  & european_subset
    subset_younger <- !is.na(pheno[, predictor]) & pheno$agedeath <= median_age & pheno$CTE == '1'  & european_subset
    # ---- older:
    fit <- with(imp, polr(as.formula(model1.form), na.action="na.omit", Hess = T, subset = subset_older))
    pool.fit <- pool(fit)
    result <- summary(pool.fit)
    assoc_subject <- !(is.na(complete(imp,1)[,.pheno])) & subset_older
    n <- sum(assoc_subject)
    effect <- round(result$estimate[result$term == predictor], 3)
    se <- round(result$std.error[result$term == predictor], 3)
    p  <- result$p.value[result$term == predictor]
    freq <- round(sum(pheno[,predictor][assoc_subject])/(2*n), 3)
    output.row <- c('Cases', 'Older', '1', .pheno, predictor, freq, n, effect, se, p)
    output.row <- unname(output.row)
    model.age.data <- rbind(model.age.data, output.row)
    # ---- younger:
    fit <- with(imp, polr(as.formula(model1.form), na.action="na.omit", Hess = T, subset = subset_younger))
    pool.fit <- pool(fit)
    result <- summary(pool.fit)
    assoc_subject <- !(is.na(complete(imp,1)[,.pheno])) & subset_younger
    n <- sum(assoc_subject)
    effect <- round(result$estimate[result$term == predictor], 3)
    se <- round(result$std.error[result$term == predictor], 3)
    p  <- result$p.value[result$term == predictor]
    freq <- round(sum(pheno[,predictor][assoc_subject])/(2*n), 3)
    output.row <- c('Cases', 'Younger', '1', .pheno, predictor, freq, n, effect, se, p)
    output.row <- unname(output.row)
    model.age.data <- rbind(model.age.data, output.row)
    
    # Fit Two All:
    subset_older <- !is.na(pheno[, predictor]) & pheno$agedeath > median_age & !is.na(pheno$footyrs)  & european_subset
    subset_younger <- !is.na(pheno[, predictor]) & pheno$agedeath <= median_age & !is.na(pheno$footyrs)  & european_subset
    # ---- older:
    fit <- with(imp, polr(as.formula(model2.form), na.action="na.omit", Hess = T, subset = subset_older))
    pool.fit <- pool(fit)
    result <- summary(pool.fit)
    assoc_subject <- !(is.na(complete(imp,1)[,.pheno])) & subset_older
    n <- sum(assoc_subject)
    effect <- round(result$estimate[result$term == predictor], 3)
    se <- round(result$std.error[result$term == predictor], 3)
    p  <- result$p.value[result$term == predictor]
    freq <- round(sum(pheno[,predictor][assoc_subject])/(2*n), 3)
    output.row <- c('All', 'Older', '2', .pheno, predictor, freq, n, effect, se, p)
    output.row <- unname(output.row)
    model.age.data <- rbind(model.age.data, output.row)
    # ---- younger:
    fit <- with(imp, polr(as.formula(model2.form), na.action="na.omit", Hess = T, subset = subset_younger))
    pool.fit <- pool(fit)
    result <- summary(pool.fit)
    assoc_subject <- !(is.na(complete(imp,1)[,.pheno])) & subset_younger
    n <- sum(assoc_subject)
    effect <- round(result$estimate[result$term == predictor], 3)
    se <- round(result$std.error[result$term == predictor], 3)
    p  <- result$p.value[result$term == predictor]
    freq <- round(sum(pheno[,predictor][assoc_subject])/(2*n), 3)
    output.row <- c('All', 'Younger', '2', .pheno, predictor, freq, n, effect, se, p)
    output.row <- unname(output.row)
    model.age.data <- rbind(model.age.data, output.row)
    # Fit Two Cases:
    subset_older <- !is.na(pheno[, predictor]) & pheno$agedeath > median_age & pheno$CTE == '1' & !is.na(pheno$footyrs) & european_subset
    subset_younger <- !is.na(pheno[, predictor]) & pheno$agedeath <= median_age & pheno$CTE == '1' & !is.na(pheno$footyrs) & european_subset
    # ---- older:
    fit <- with(imp, polr(as.formula(model2.form), na.action="na.omit", Hess = T, subset = subset_older))
    pool.fit <- pool(fit)
    result <- summary(pool.fit)
    assoc_subject <- !(is.na(complete(imp,1)[,.pheno])) & subset_older
    n <- sum(assoc_subject)
    effect <- round(result$estimate[result$term == predictor], 3)
    se <- round(result$std.error[result$term == predictor], 3)
    p  <- result$p.value[result$term == predictor]
    freq <- round(sum(pheno[,predictor][assoc_subject])/(2*n), 3)
    output.row <- c('Cases', 'Older', '2', .pheno, predictor, freq, n, effect, se, p)
    output.row <- unname(output.row)
    model.age.data <- rbind(model.age.data, output.row)
    # ---- younger:
    fit <- with(imp, polr(as.formula(model2.form), na.action="na.omit", Hess = T, subset = subset_younger))
    pool.fit <- pool(fit)
    result <- summary(pool.fit)
    assoc_subject <- !(is.na(complete(imp,1)[,.pheno])) & subset_younger
    n <- sum(assoc_subject)
    effect <- round(result$estimate[result$term == predictor], 3)
    se <- round(result$std.error[result$term == predictor], 3)
    p  <- result$p.value[result$term == predictor]
    freq <- round(sum(pheno[,predictor][assoc_subject])/(2*n), 3)
    output.row <- c('Cases', 'Younger', '2', .pheno, predictor, freq, n, effect, se, p)
    output.row <- unname(output.row)
    model.age.data <- rbind(model.age.data, output.row)
    
    # Fit Three All:
    subset_older <- !is.na(pheno[, predictor]) & pheno$agedeath > median_age & !is.na(pheno$footyrs) & european_subset
    subset_younger <- !is.na(pheno[, predictor]) & pheno$agedeath <= median_age & !is.na(pheno$footyrs) & european_subset
    # ---- older:
    fit <- with(imp, polr(as.formula(model3.form), na.action="na.omit", Hess = T, subset = subset_older))
    pool.fit <- pool(fit)
    result <- summary(pool.fit)
    assoc_subject <- !(is.na(complete(imp,1)[,.pheno])) & subset_older
    n <- sum(assoc_subject)
    effect <- round(result$estimate[result$term == paste0(predictor, ":footyrs")], 3)
    se <- round(result$std.error[result$term == paste0(predictor, ":footyrs")], 3)
    p  <- result$p.value[result$term == paste0(predictor, ":footyrs")]
    freq <- round(sum(pheno[,predictor][assoc_subject])/(2*n), 3)
    output.row <- c('All', 'Older', '3', .pheno, paste0(predictor, "*footyrs"), freq, n, effect, se, p)
    output.row <- unname(output.row)
    model.age.data <- rbind(model.age.data, output.row)
    # ---- younger:
    fit <- with(imp, polr(as.formula(model3.form), na.action="na.omit", Hess = T, subset = subset_younger))
    pool.fit <- pool(fit)
    result <- summary(pool.fit)
    assoc_subject <- !(is.na(complete(imp,1)[,.pheno])) & subset_younger
    n <- sum(assoc_subject)
    effect <- round(result$estimate[result$term == paste0(predictor, ":footyrs")], 3)
    se <- round(result$std.error[result$term == paste0(predictor, ":footyrs")], 3)
    p  <- result$p.value[result$term == paste0(predictor, ":footyrs")]
    freq <- round(sum(pheno[,predictor][assoc_subject])/(2*n), 3)
    output.row <- c('All', 'Younger', '3', .pheno, paste0(predictor, "*footyrs"), freq, n, effect, se, p)
    output.row <- unname(output.row)
    model.age.data <- rbind(model.age.data, output.row)
    # Fit Three Cases:
    subset_older <- !is.na(pheno[, predictor]) & pheno$agedeath > median_age & pheno$CTE == '1' & !is.na(pheno$footyrs) & european_subset
    subset_younger <- !is.na(pheno[, predictor]) & pheno$agedeath <= median_age & pheno$CTE == '1' & !is.na(pheno$footyrs) & european_subset
    # ---- older:
    fit <- with(imp, polr(as.formula(model3.form), na.action="na.omit", Hess = T, subset = subset_older))
    pool.fit <- pool(fit)
    result <- summary(pool.fit)
    assoc_subject <- !(is.na(complete(imp,1)[,.pheno])) & subset_older
    n <- sum(assoc_subject)
    effect <- round(result$estimate[result$term == paste0(predictor, ":footyrs")], 3)
    se <- round(result$std.error[result$term == paste0(predictor, ":footyrs")], 3)
    p  <- result$p.value[result$term == paste0(predictor, ":footyrs")]
    freq <- round(sum(pheno[,predictor][assoc_subject])/(2*n), 3)
    output.row <- c('Cases', 'Older', '3', .pheno, paste0(predictor, "*footyrs"), freq, n, effect, se, p)
    output.row <- unname(output.row)
    model.age.data <- rbind(model.age.data, output.row)
    # ---- younger:
    fit <- with(imp, polr(as.formula(model3.form), na.action="na.omit", Hess = T, subset = subset_younger))
    pool.fit <- pool(fit)
    result <- summary(pool.fit)
    assoc_subject <- !(is.na(complete(imp,1)[,.pheno])) & subset_younger
    n <- sum(assoc_subject)
    effect <- round(result$estimate[result$term == paste0(predictor, ":footyrs")], 3)
    se <- round(result$std.error[result$term == paste0(predictor, ":footyrs")], 3)
    p  <- result$p.value[result$term == paste0(predictor, ":footyrs")]
    freq <- round(sum(pheno[,predictor][assoc_subject])/(2*n), 3)
    output.row <- c('Cases', 'Younger', '3', .pheno, paste0(predictor, "*footyrs"), freq, n, effect, se, p)
    output.row <- unname(output.row)
    model.age.data <- rbind(model.age.data, output.row)
    
    # Fit Four All:
    subset <- !is.na(pheno[, predictor])
    fit <- with(imp, polr(as.formula(model4.form), na.action="na.omit", Hess = T, subset = subset))
    pool.fit <- pool(fit)
    result <- summary(pool.fit)
    assoc_subject <- !(is.na(complete(imp,1)[,.pheno])) & subset
    n <- sum(assoc_subject)
    effect <- round(result$estimate[result$term == paste0(predictor, ":stratum")], 3)
    se <- round(result$std.error[result$term == paste0(predictor, ":stratum")], 3)
    p  <- result$p.value[result$term == paste0(predictor, ":stratum")]
    freq <- round(sum(pheno[,predictor][assoc_subject])/(2*n), 3)
    output.row <- c('All', 'NA', '4', .pheno, paste0(predictor, "*stratum"), freq, n, effect, se, p)
    output.row <- unname(output.row)
    model.age.data <- rbind(model.age.data, output.row)
    # Fit Four Cases:
    subset <- !is.na(pheno[, predictor]) & pheno$CTE == '1'
    fit <- with(imp, polr(as.formula(model4.form), na.action="na.omit", Hess = T, subset = subset))
    pool.fit <- pool(fit)
    result <- summary(pool.fit)
    assoc_subject <- !(is.na(complete(imp,1)[,.pheno])) & subset
    n <- sum(assoc_subject)
    effect <- round(result$estimate[result$term == paste0(predictor, ":stratum")], 3)
    se <- round(result$std.error[result$term == paste0(predictor, ":stratum")], 3)
    p  <- result$p.value[result$term == paste0(predictor, ":stratum")]
    freq <- round(sum(pheno[,predictor][assoc_subject])/(2*n), 3)
    output.row <- c('Cases', 'NA', '4', .pheno, paste0(predictor, "*stratum"), freq, n, effect, se, p)
    output.row <- unname(output.row)
    model.age.data <- rbind(model.age.data, output.row)
  }
  # Calc FDR
  model.age.data <- model.age.data[2:nrow(model.age.data), ]
  model.age.data$FDR <- NA
  model.age.data$FDR[model.age.data$Model == '1'] <- p.adjust(model.age.data$P.value[model.age.data$Model == '1'], method = 'fdr')
  model.age.data$FDR[model.age.data$Model == '2'] <- p.adjust(model.age.data$P.value[model.age.data$Model == '2'], method = 'fdr')
  model.age.data$FDR[model.age.data$Model == '3'] <- p.adjust(model.age.data$P.value[model.age.data$Model == '3'], method = 'fdr')
  model.age.data$FDR[model.age.data$Model == '4'] <- p.adjust(model.age.data$P.value[model.age.data$Model == '4'], method = 'fdr')
  # save result
  outputname <- 'results/Model.Age_Stratified.XXXX.assoc.txt'
  outputname <- gsub('XXXX', predictor, outputname)
  write.table(model.age.data, outputname, col.names = T, row.names = F, sep = '\t', quote = F)
}
