# mice imputation including AT8_positive_tissue and run the association test
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
at8pt <- haven::read_sav('/restricted/projectnb/cte/Challenge_Project_2022/January 2023 AT8 Update/APOEAnalysiswithupdatedAT8data.sav')

pheno$AT8PositiveTissue <- at8pt$AT8PositiveTissue[match(pheno$subjid, at8pt$subjid)]
pheno$AT8PositiveTissue.log <- log(pheno$AT8PositiveTissue)
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

# set variables
pheno_list_categorical <- c('amygdala','CA1','CA2','CA4','entorhinal','inforbfront',
                            'infpar','loccoer','midfront','suptemp','CTEstage','subnigra')

pheno_list_binary      <- c('CTE', 'DementiaHx')

pheno_list_continuous  <- c('AT8PositiveTissue.log')

pheno_other_info       <- c('pin', 'subjid', 'agedeath', 'race', 'footyrs')

# keep the useful phenotypes 
pheno <- pheno[, c(pheno_other_info, pheno_list_binary, pheno_list_categorical, pheno_list_continuous)]

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
batch2.subhap$ID[batch2.subhap$ID == 'k-0487'] <- 'K-0487'
batch2.subhap$ID[batch2.subhap$ID == 'k-0554'] <- 'K-0554'
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

unite2.pcs$V1[unite2.pcs$V1 == 'k-0487'] <- 'K-0487'
unite2.pcs$V1[unite2.pcs$V1 == 'k-0554'] <- 'K-0554'
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

# convert H1B1G1 to dominant model
pheno$H1B1G1 <- ifelse(pheno$H1B1G1 == 2, 1, pheno$H1B1G1)
# Imputation 
# =============================================================================================
# initiation
ini <- mice(pheno, maxit=0, print=F)
pred <- ini$predictorMatrix
pred[!(rownames(pred) %in% c(pheno_list_categorical, pheno_list_continuous)),] <- 0
pred['CTEstage',] <- 0
pred[,!(colnames(pred) %in% c(pheno_list_categorical,'agedeath', pheno_list_continuous))] <- 0
pred[,'CTEstage'] <- 0

# try again
ini <- mice(pheno, maxit=0, print=F, pred = pred)
method <- ini$method
method[pheno_list_categorical] <- "polr"
method[pheno_list_continuous]  <- "norm"
method[!names(method) %in% c(pheno_list_categorical, pheno_list_continuous)] <- ""
method['CTEstage'] <- ""

# run the imputation
imp <-  mice(pheno, pred = pred, m = m, method = method, 
             maxit = maxit, seed = seed, print = T)
# =============================================================================================

# save the imputation result
save.image('imputation.RData')
