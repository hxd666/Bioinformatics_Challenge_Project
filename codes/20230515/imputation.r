# mice imputation including AT8_positive_tissue and run the association test
rm(list = ls())
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
pheno$subjid[pheno$subjid == 'k-0487'] <- "K-0487"
pheno$subjid[pheno$subjid == 'k-0554'] <- 'K-0554'
# read the at8 positive tissue data
at8pt <- haven::read_sav('/restricted/projectnb/cte/Challenge_Project_2022/data/pheno/update_05152023/APOEAnalysisupdatedAT8data.sav')

pheno$AT8PositiveTissue <- at8pt$AT8PositiveTissue[match(pheno$subjid, at8pt$subjid)]
pheno$AT8PositiveTissue.log <- log(pheno$AT8PositiveTissue)
pheno$HipCA23 <- at8pt$HipCA23[match(pheno$subjid, at8pt$subjid)]
pheno$HipCA23.log <- log(pheno$HipCA23+1e-05)
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

# check log AT8 positive tissue (missing 81)
table(pheno$AT8PositiveTissue.log, useNA = 'always')

# check quantitative CA2/CA3 (missing 136)
table(pheno$HipCA23.log, useNA = 'always')

# set variables
pheno_list_categorical <- c('amygdala','CA1','CA2','CA4','entorhinal','inforbfront',
                            'infpar','loccoer','midfront','suptemp','CTEstage','subnigra')

pheno_list_binary      <- c('CTE', 'DementiaHx')

pheno_list_continuous  <- c('AT8PositiveTissue.log', 'HipCA23.log')

pheno_other_info       <- c('pin', 'subjid', 'agedeath', 'race', 'footyrs')

# keep the useful phenotypes 
pheno <- pheno[, c(pheno_other_info, pheno_list_binary, pheno_list_categorical, pheno_list_continuous)]

# add dosage information
# batch1.subhap <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/imputation/1st_batch/subhaplotypes/1st_batch.topmed.mapt.subhaplotypes.txt',
#                             header = T, as.is = T)
# batch2.subhap <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/imputation/2nd_batch/subhaplotypes/2nd_batch.topmed.mapt.subhaplotypes.txt',
#                             header = T, as.is = T)
# batch1.subhap$TOPMed_ID <- sapply(batch1.subhap$TOPMed_ID, function(x){strsplit(x,'_')[[1]][1]})
# batch1.subhap$ID <- batch1.subhap$TOPMed_ID
# batch1.subhap <- batch1.subhap[batch1.subhap$ID %in% pheno$pin, ]
# batch1.subhap$Batch <- 'I'
# batch1.subhap$ID <- pheno$subjid[match(batch1.subhap$ID, pheno$pin)]


# batch2.subhap$ID <- batch2.subhap$TOPMed_ID
# batch2.subhap$ID[batch2.subhap$ID == 'k-0487'] <- 'K-0487'
# batch2.subhap$ID[batch2.subhap$ID == 'k-0554'] <- 'K-0554'
# batch2.subhap <- batch2.subhap[batch2.subhap$ID %in% pheno$subjid, ]
# batch2.subhap$Batch <- 'II'
# 
# subhap <- rbind(batch1.subhap, batch2.subhap)
# subhap <- subhap[match(pheno$subjid, subhap$ID), ]
# pheno <- cbind(pheno, subhap)
subhap <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/imputation/neural_net/Prediction_data_updated_PredResult_LATEST.txt',
                     header = T, sep = ',', as.is = T)
subhap$hap1.h1b1g1 <- ifelse(subhap$hap1 == 'H1B1G1', 1, 0)
subhap$hap2.h1b1g1 <- ifelse(subhap$hap2 == 'H1B1G1', 1, 0)
subhap$H1B1G1 <- subhap$hap1.h1b1g1 + subhap$hap2.h1b1g1
table(subhap$H1B1G1)
subhap <- subhap[, c("ID", "H1B1G1")]
subhap <- subhap[order(subhap$ID), ]

batch1.id <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/_be562/data/pred1.fam', as.is = T)
batch1.id <- data.frame('ID' = sapply(batch1.id$V2, function(x){strsplit(x,split="_")[[1]][1]}), stringsAsFactors = F)
batch1.id$Batch <- 'I'
batch1.id$H1B1G1 <- subhap$H1B1G1[1:nrow(batch1.id)]
batch2.id <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/_be562/data/pred2.fam', as.is = T)
batch2.id <- data.frame('ID' = batch2.id$V2, stringsAsFactors = F)
batch2.id$Batch <- 'II'
batch2.id$ID[batch2.id$ID == 'k-0487'] <- "K-0487"
batch2.id$ID[batch2.id$ID == 'k-0554'] <- 'K-0554'
batch2.id$H1B1G1 <- subhap$H1B1G1[(nrow(batch1.id)+1):nrow(subhap)]

batch1.id <- batch1.id[batch1.id$ID %in% pheno$pin, ]
batch1.id <- batch1.id[match(pheno$pin, batch1.id$ID), ]
batch1.id <- na.omit(batch1.id)
batch2.id <- batch2.id[batch2.id$ID %in% pheno$subjid, ]
batch2.id <- batch2.id[match(pheno$subjid, batch2.id$ID), ]
batch2.id <- na.omit(batch2.id)
# check ID
all(batch1.id$ID == pheno$pin[1:nrow(batch1.id)])
all(batch2.id$ID == pheno$subjid[(nrow(batch1.id)+1):nrow(pheno)])
subhap <- rbind(batch1.id, batch2.id)

pheno <- cbind(pheno, subhap)

batch1.id <- subhap[subhap$Batch == 'I', ]
batch2.id <- subhap[subhap$Batch == 'II', ]
# Add PCs
unite1.pcs <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/UNITE_batch1/PCA/UNITE_batch1.pca.txt', header = F, as.is = T)
unite2.pcs <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/UNITE_batch2/PCA/UNITE_batch2.pca.txt', header = F, as.is = T)

unite1.pcs <- unite1.pcs[unite1.pcs$V1 %in% batch1.id$ID, ]
unite1.pcs$V1 <- batch1.id$ID[match(unite1.pcs$V1, batch1.id$ID)]
colnames(unite1.pcs) <- c('ID', sapply(1:20, function(x){paste0('PC',x)}))

unite2.pcs$V1[unite2.pcs$V1 == 'k-0487'] <- 'K-0487'
unite2.pcs$V1[unite2.pcs$V1 == 'k-0554'] <- 'K-0554'
unite2.pcs <- unite2.pcs[unite2.pcs$V1 %in% batch2.id$ID, ]
colnames(unite2.pcs) <- c('ID', sapply(1:20, function(x){paste0('PC',x)}))

pcs <- rbind(unite1.pcs, unite2.pcs)
pcs <- pcs[match(pheno$ID, pcs$ID), ]
pheno <- merge(pheno, pcs, by = 'ID')

# add stratum
median_age <- median(pheno$agedeath)
pheno$stratum <- ifelse(pheno$agedeath <= median_age, 0, 1)

# add rs6910507
rs_1st <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20230515/1st_batch.rs6910507.dose', as.is = T)
rs_1st$V1 <- NULL; rs_1st$V2 <- NULL
rs_1st_list <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20230515/1st_sample', as.is = T)
rs_1st <- cbind(rs_1st_list, t(rs_1st))
colnames(rs_1st) <- c('ID', 'rs6910507')
rs_1st$ID <- sapply(rs_1st$ID, function(x){strsplit(x,'_')[[1]][1]})

rs_2nd <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20230515/2nd_batch.rs6910507.dose', as.is = T)
rs_2nd$V1 <- NULL; rs_2nd$V2 <- NULL
rs_2nd_list <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20230515//2nd_sample', as.is = T)
rs_2nd <- cbind(rs_2nd_list, t(rs_2nd))
colnames(rs_2nd) <- c('ID', 'rs6910507')

rs <- rbind(rs_1st, rs_2nd)
rs$ID[rs$ID == 'k-0487'] <- 'K-0487'
rs$ID[rs$ID == 'k-0554'] <- 'K-0554'
sum(pheno$ID %in% rs$ID) == nrow(pheno)
pheno$rs6910507 <- rs$rs6910507[match(pheno$ID, rs$ID)]

# convert the categorical phenotypes to factor datatype
for (.pheno in pheno_list_categorical){
  pheno[, .pheno] <- as.factor(pheno[, .pheno])
}
# Imputation 
# =============================================================================================
# initiation
ini <- mice(pheno, maxit=0, print=F)
pred <- ini$predictorMatrix
pred[!(rownames(pred) %in% c(pheno_list_categorical, pheno_list_continuous)),] <- 0
pred['CTEstage',] <- 0
pred[,!(colnames(pred) %in% c(pheno_list_categorical,'agedeath', pheno_list_continuous))] <- 0
pred[,'CTEstage'] <- 0
pred[,'HipCA23.log'] <- 0


# try again
ini <- mice(pheno, maxit=0, print=F, pred = pred)
method <- ini$method
method[pheno_list_categorical] <- "polr"
method[pheno_list_continuous]  <- "norm"
method[!names(method) %in% c(pheno_list_categorical, pheno_list_continuous)] <- ""
method['CTEstage'] <- ""
#method['HipCA23.log'] <- ""
# run the imputation
imp <-  mice(pheno, pred = pred, m = m, method = method, 
             maxit = maxit, seed = seed, print = T)
# =============================================================================================

# save the imputation result
save.image('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20230515/imputation.RData')

