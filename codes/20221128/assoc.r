# Do the association test
require(mice)
require(MASS)
require(rms)
# =========================
# Mice Imputation Parameter
# =========================
## m: Number of data set we want
## seed: Set number for replication of simulation
## print: additional info you can have printed
## maxit: Number of iterations in the MCMC
## donors: Number of nearest neighbors to define the PMM
## method: Defines the type of variable to impute
m <- 10
seed <- 12012020
maxit <- 20
donors <- 5

# Curate the phenotype file:
pheno <- read.csv('/restricted/projectnb/cte/Challenge_Project_2022/data/pheno/update_11182022/APOEAnalysis_UNITE_wo_unmatched.csv',
                  as.is = T)
# keep the sample with age death info
pheno <- pheno[complete.cases(pheno$agedeath), ]
# keep those age >= 20
pheno <- pheno[pheno$agedeath >= 20, ]

# check CTE (no missing)
table(pheno$CTE, useNA = 'always')
pheno$CTE[pheno$CTE=='Yes'] <- "1"
pheno$CTE <- as.numeric(pheno$CTE)
table(pheno$CTE, useNA = 'always')

# check CTE stage (missing one)
table(pheno$CTEstage, useNA = 'always')
pheno$CTEstage[pheno$CTEstage %in% 'Stage IV'] <- '4'
pheno$CTEstage <- as.numeric(pheno$CTEstage)
table(pheno$CTEstage, useNA = 'always')

# check Amygdala (missing 10)
table(pheno$amygdala, useNA = 'always')
pheno$amygdala <- as.numeric(pheno$amygdala)

# check CA1 (missing 10)
table(pheno$CA1, useNA = 'always')
pheno$CA1 <- as.numeric(pheno$CA1)

# check CA2 (missing 31)
table(pheno$CA2, useNA = 'always')
pheno$CA2 <- as.numeric(pheno$CA2)

# check CA4 (missing 16)
table(pheno$CA4, useNA = 'always')
pheno$CA4 <- as.numeric(pheno$CA4)

# check Entorhinal (missing 8)
table(pheno$entorhinal, useNA = 'always')
pheno$entorhinal <- as.numeric(pheno$entorhinal)

# check InfOrbFront (missing 10)
table(pheno$inforbfront, useNA = 'always')
pheno$inforbfront <- as.numeric(pheno$inforbfront)

# check InfPar (missing 8)
table(pheno$infpar, useNA = 'always')
pheno$infpar <- as.numeric(pheno$infpar)

# check LocCoer (missing 32)
table(pheno$loccoer, useNA = 'always')
pheno$loccoer <- as.numeric(pheno$loccoer)

# check MidFront (missing 1)
table(pheno$midfront, useNA = 'always')
pheno$midfront <- as.numeric(pheno$midfront)

# check SupTemp (missing 8)
table(pheno$suptemp, useNA = 'always')
pheno$suptemp <- as.numeric(pheno$suptemp)

# check SubNigra (missing 7)
table(pheno$subnigra, useNA = 'always')
pheno$subnigra <- as.numeric(pheno$subnigra)

# check race
table(pheno$race, useNA = 'always')
pheno$race[pheno$race %in% 'White'] <- '1'
table(pheno$race, useNA = 'always')

# check dementia (missing 32)
table(pheno$DementiaHx, useNA = 'always')
pheno$DementiaHx[pheno$DementiaHx %in% 'Yes'] <- '1'
pheno$DementiaHx <- as.numeric(pheno$DementiaHx)

# set variables
pheno_list_categorical <- c('amygdala','CA1','CA2','CA4','entorhinal','inforbfront',
                            'infpar','loccoer','midfront','suptemp','CTEstage','subnigra')

pheno_list_binary      <- c('CTE', 'DementiaHx')

pheno_other_info       <- c('pin', 'subjid', 'agedeath', 'race', 'footyrs')

# keep the useful phenotypes 
pheno <- pheno[, c(pheno_other_info, pheno_list_binary, pheno_list_categorical)]

# add dosage information
batch1.subhap <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/imputation/1st_batch/subhaplotypes/1st_batch.topmed.mapt.subhaplotypes.txt',
                            header = T, as.is = T)
batch2.subhap <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/imputation/2nd_batch/subhaplotypes/2nd_batch.topmed.mapt.subhaplotypes.txt',
                            header = T, as.is = T)
batch1.subhap$TOPMed_ID <- sapply(batch1.subhap$TOPMed_ID, function(x){strsplit(x,'_')[[1]][1]})
batch1.subhap$ID <- batch1.subhap$TOPMed_ID
batch1.subhap <- batch1.subhap[batch1.subhap$ID %in% pheno$pin, ]
batch1.subhap$Batch <- 'I'
batch1.subhap$ID <- pheno$subjid[match(batch1.subhap$ID, pheno$pin)]


batch2.subhap$ID <- batch2.subhap$TOPMed_ID
batch2.subhap <- batch2.subhap[batch2.subhap$ID %in% pheno$subjid, ]
batch2.subhap$Batch <- 'II'

subhap <- rbind(batch1.subhap, batch2.subhap)
subhap <- subhap[match(pheno$subjid, subhap$ID), ]
pheno <- cbind(pheno, subhap)

# Add PCs
unite1.pcs <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/UNITE_batch1/PCA/UNITE_batch1.pca.txt', header = F, as.is = T)
unite2.pcs <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/UNITE_batch2/PCA/UNITE_batch2.pca.txt', header = F, as.is = T)

unite1.pcs <- unite1.pcs[unite1.pcs$V1 %in% batch1.subhap$TOPMed_ID, ]
unite1.pcs$V1 <- batch1.subhap$ID[match(unite1.pcs$V1, batch1.subhap$TOPMed_ID)]
colnames(unite1.pcs) <- c('ID', sapply(1:20, function(x){paste0('PC',x)}))

unite2.pcs <- unite2.pcs[unite2.pcs$V1 %in% batch2.subhap$ID, ]
colnames(unite2.pcs) <- c('ID', sapply(1:20, function(x){paste0('PC',x)}))

pcs <- rbind(unite1.pcs, unite2.pcs)
pcs <- pcs[match(pheno$ID, pcs$ID), ]
pheno <- merge(pheno, pcs, by = 'ID')

# add stratum
median_age <- median(pheno$agedeath)
pheno$stratum <- ifelse(pheno$agedeath <= median_age, 0, 1)

# convert the categorical phenotypes to factor datatype
for (.pheno in pheno_list_categorical){
  pheno[, .pheno] <- as.factor(pheno[, .pheno])
}
# Imputation 
# =============================================================================================
# initiation
ini <- mice(pheno, maxit=0, print=F)
pred <- ini$predictorMatrix
pred[!(rownames(pred) %in% c(pheno_list_categorical)),] <- 0
pred['CTEstage',] <- 0
pred[,!(colnames(pred) %in% c(pheno_list_categorical,'agedeath'))] <- 0
pred[,'CTEstage'] <- 0

# try again
ini <- mice(pheno, maxit=0, print=F, pred = pred)
method <- ini$method
method[pheno_list_categorical] <- "polyreg"
method[!names(method) %in% pheno_list_categorical] <- ""
method['CTEstage'] <- ""

# run the imputation
imp <-  mice(pheno, pred = pred, m = m, method = method, 
             maxit = maxit, seed = seed, print = T)
# =============================================================================================

# save the imputation result
save.image('imputatoin.RData')


# Check the Principal Components
pcs.formula <- sapply(1:10, function(x){paste0('PC',x)})
pcs.formula <- paste(pcs.formula, collapse = ' + ')

model0 <- function(pheno1){
  form <- paste(pheno1, " ~  ", pcs.formula, sep = '')
  fit <- with(imp, polr(as.formula(form), na.action="na.omit", Hess = T))
  pool.fit <- pool(fit)
  result <- summary(pool.fit)
  n <- sum(!is.na(complete(imp,1)[,pheno1]))
  p.values <- result$p.value[1:10]
  output.row <- c( pheno1, n , p.values)
  output.row <- unname(output.row)
  return(output.row)
}
model0.data <- data.frame('Phenotype' = 0, 'N' = 0,
                          'PC1' = 0, 'PC2' = 0, 'PC3' = 0, 'PC4' = 0,
                          'PC5' = 0, 'PC6' = 0, 'PC7' = 0, 'PC8' = 0,
                          'PC9' = 0, 'PC10' = 0, stringsAsFactors = F)

for (.pheno in pheno_list_categorical){
  print(.pheno)
  model0.data <- rbind(model0.data, model0(.pheno))
}
model0.data <- model0.data[2:nrow(model0.data), ]
model0.data[, 3:12] <- apply(model0.data[, 3:12], 2, function(x){p.adjust(x, method = 'fdr')})
write.table(model0.data, 'model0.euro.txt', col.names = T, row.names = F, sep = '\t', quote = F)




# model 1 Phenotype ~ H1_dosage + Age + PCs 
euro.batch1 <- pheno$race %in% 1 & pheno$stratum == 0
euro.batch2 <- pheno$race %in% 1 & pheno$stratum == 1
model1 <- function(pheno1, batch){
  form <- paste(pheno1, " ~ H1_dosage + agedeath + ", pcs.formula, sep = '')
  fit <- with(imp, polr(as.formula(form), na.action="na.omit", Hess = T, subset = batch))
  pool.fit <- pool(fit)
  result <- summary(pool.fit)
  if (all(batch == euro.batch1)){
    batch.n = '0'
    } else {
      batch.n = '1'
    }
  n <- sum(!is.na(complete(imp,1)[batch,pheno1]))
  effect <- result$estimate[result$term == 'H1_dosage']
  se <- result$std.error[result$term == 'H1_dosage']
  p  <- result$p.value[result$term == 'H1_dosage']
  output.row <- c(batch.n, pheno1, 'H1_dosage', n, effect, se, p)
  output.row <- unname(output.row)
  return(output.row)
}
model1.data <- data.frame('Stratum' = 0,'Phenotype' = 0, 'Predictor'= 0, 'N' = 0,
                          'Effect' = 0, 'SE' = 0, 'P-value' = 0, stringsAsFactors = F)
for (.pheno in pheno_list_categorical){
  print(.pheno)
  model1.data <- rbind(model1.data, model1(.pheno, euro.batch1))
  model1.data <- rbind(model1.data, model1(.pheno, euro.batch2))
}
model1.data <- model1.data[2:nrow(model1.data), ]
model1.data$FDR <- p.adjust(model1.data$P.value, method = 'fdr')
write.table(model1.data, 'model1.euro.txt', col.names = T, row.names = F, sep = '\t', quote = F)


# model 2 Phenotype ~ H1_dosage + Age + PCs cases only
euro.batch1.cases <- pheno$race %in% 1 & pheno$CTE == 1 & pheno$stratum == 0
euro.batch2.cases <- pheno$race %in% 1 & pheno$CTE == 1 & pheno$stratum == 1
model2 <- function(pheno1, batch){
  form <- paste(pheno1, " ~ H1_dosage + agedeath + ", pcs.formula, sep = '')
  fit <- with(imp, polr(as.formula(form), na.action="na.omit", Hess = T, subset = batch))
  pool.fit <- pool(fit)
  result <- summary(pool.fit)
  if (all(batch == euro.batch1.cases)){
    batch.n = '0'
    } else {
      batch.n = '1'
    }
  n <- sum(!is.na(complete(imp,1)[batch,pheno1]))
  effect <- result$estimate[result$term == 'H1_dosage']
  se <- result$std.error[result$term == 'H1_dosage']
  p  <- result$p.value[result$term == 'H1_dosage']
  output.row <- c(batch.n, pheno1, 'H1_dosage', n, effect, se, p)
  output.row <- unname(output.row)
  return(output.row)
}

model2.data <- data.frame('Stratum' = 0,'Phenotype' = 0, 'Predictor'= 0, 'N' = 0,
                          'Effect' = 0, 'SE' = 0, 'P-value' = 0, stringsAsFactors = F)
for (.pheno in pheno_list_categorical){
  print(.pheno)
  model2.data <- rbind(model2.data, model2(.pheno, euro.batch1.cases))
  model2.data <- rbind(model2.data, model2(.pheno, euro.batch2.cases))
}
model2.data <- model2.data[2:nrow(model2.data), ]
model2.data$FDR <- p.adjust(model2.data$P.value, method = 'fdr')
write.table(model2.data, 'model2.euro.txt', col.names = T, row.names = F, sep = '\t', quote = F)


# model 3 Phenotype ~ H1_dosage + Age + PCs cases only + Duration in football players only 
euro.batch1.football <- pheno$race %in% 1 & !is.na(pheno$footyrs) & pheno$stratum == 0
euro.batch2.football <- pheno$race %in% 1 & !is.na(pheno$footyrs) & pheno$stratum == 1
model3 <- function(pheno1, batch){
  form <- paste(pheno1, " ~ H1_dosage + agedeath + footyrs + ", pcs.formula, sep = '')
  fit <- with(imp, polr(as.formula(form), na.action="na.omit", Hess = T, subset = batch))
  pool.fit <- pool(fit)
  result <- summary(pool.fit)
  if (all(batch == euro.batch1.football)){
    batch.n = '0'
  } else {
    batch.n = '1'
  }
  n <- sum(!is.na(complete(imp,1)[batch,pheno1]))
  effect <- result$estimate[result$term == 'H1_dosage']
  se <- result$std.error[result$term == 'H1_dosage']
  p  <- result$p.value[result$term == 'H1_dosage']
  output.row <- c(batch.n, pheno1, 'H1_dosage', n, effect, se, p)
  output.row <- unname(output.row)
  return(output.row)
}

model3.data <- data.frame('Stratum' = 0,'Phenotype' = 0, 'Predictor'= 0, 'N' = 0,
                          'Effect' = 0, 'SE' = 0, 'P-value' = 0, stringsAsFactors = F)
for (.pheno in pheno_list_categorical){
  print(.pheno)
  model3.data <- rbind(model3.data, model3(.pheno, euro.batch1.football))
  model3.data <- rbind(model3.data, model3(.pheno, euro.batch2.football))
}
model3.data <- model3.data[2:nrow(model3.data), ]
model3.data$FDR <- p.adjust(model3.data$P.value, method = 'fdr')
write.table(model3.data, 'model3.euro.txt', col.names = T, row.names = F, sep = '\t', quote = F)

# model 4 Phenotype ~ H1_dosage + Age + PCs cases only + Duration + H1_dosage*Duration in football players only

model4 <- function(pheno1, batch){
  form <- paste(pheno1, " ~ H1_dosage + agedeath + footyrs + H1_dosage*footyrs+ ", pcs.formula, sep = '')
  fit <- with(imp, polr(as.formula(form), na.action="na.omit", Hess = T, subset = batch))
  pool.fit <- pool(fit)
  result <- summary(pool.fit)
  if (all(batch == euro.batch1.football)){
    batch.n = '0'
  } else {
    batch.n = '1'
  }
  n <- sum(!is.na(complete(imp,1)[batch,pheno1]))
  effect <- result$estimate[result$term == 'H1_dosage:footyrs']
  se <- result$std.error[result$term == 'H1_dosage:footyrs']
  p  <- result$p.value[result$term == 'H1_dosage:footyrs']
  output.row <- c(batch.n, pheno1, 'H1_dosage:footyrs', n, effect, se, p)
  output.row <- unname(output.row)
  return(output.row)
}

model4.data <- data.frame('Stratum' = 0,'Phenotype' = 0, 'Predictor'= 0, 'N' = 0,
                          'Effect' = 0, 'SE' = 0, 'P-value' = 0, stringsAsFactors = F)
for (.pheno in pheno_list_categorical){
  print(.pheno)
  model4.data <- rbind(model4.data, model4(.pheno, euro.batch1.football))
  model4.data <- rbind(model4.data, model4(.pheno, euro.batch2.football))
}
model4.data <- model4.data[2:nrow(model4.data), ]
model4.data$FDR <- p.adjust(model4.data$P.value, method = 'fdr')
write.table(model4.data, 'model4.euro.txt', col.names = T, row.names = F, sep = '\t', quote = F)


# model 5 Phenotype ~ H1_dosage + Stratum + H1_dosage*Stratum + PCs

model5 <- function(pheno1){
  form <- paste(pheno1, " ~ H1_dosage + stratum + H1_dosage*stratum + ", pcs.formula, sep = '')
  fit <- with(imp, polr(as.formula(form), na.action="na.omit", Hess = T))
  pool.fit <- pool(fit)
  result <- summary(pool.fit)
  n <- sum(!is.na(complete(imp,1)[,pheno1]))
  effect <- result$estimate[result$term == 'H1_dosage:stratum']
  se <- result$std.error[result$term == 'H1_dosage:stratum']
  p  <- result$p.value[result$term == 'H1_dosage:stratum']
  output.row <- c(pheno1, 'H1_dosage:stratum', n, effect, se, p)
  output.row <- unname(output.row)
  return(output.row)
}

model5.data <- data.frame('Phenotype' = 0, 'Predictor'= 0, 'N' = 0,
                          'Effect' = 0, 'SE' = 0, 'P-value' = 0, stringsAsFactors = F)
for (.pheno in pheno_list_categorical){
  print(.pheno)
  model5.data <- rbind(model5.data, model5(.pheno))
}
model5.data <- model5.data[2:nrow(model5.data), ]
model5.data$FDR <- p.adjust(model5.data$P.value, method = 'fdr')
write.table(model5.data, 'model5.euro.txt', col.names = T, row.names = F, sep = '\t', quote = F)
