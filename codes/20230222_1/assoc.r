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
load('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20230212/imputation.RData')
# Set PCs formula
#pred.list   <- read.table('predictors_list', as.is = T)
#pred.list   <- pred.list$V1


# Check the common individuals
subhap <- read.csv('/restricted/projectnb/cte/Challenge_Project_2022/imputation/neural_net/Prediction_data_updated_PredResult.txt',
                    as.is = T)

subhap <- subhap[order(subhap$ID), ]
subhap.fam1 <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/_be562/data/pred1.fam',
                          header = F, as.is = T)$V2
subhap.fam1 <- sapply(subhap.fam1, function(x){strsplit(x,split='_')[[1]][1]})
subhap.fam2 <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/_be562/data/pred2.fam',
                          header = F, as.is = T)$V2
subhap$ID <- c(subhap.fam1, subhap.fam2)
subhap.batch1 <- subhap[1:length(subhap.fam1), ]
subhap.batch2 <- subhap[(length(subhap.fam1)+1):nrow(subhap), ]
subhap.batch1 <- subhap.batch1[subhap.batch1$ID %in% pheno$pin[pheno$Batch=='I'], ]
subhap.batch2$ID[grepl('k-', subhap.batch2$ID)] <- sapply(subhap.batch2$ID[grepl('k-', subhap.batch2$ID)], function(x){gsub('k','K',x)})
subhap.batch2 <- subhap.batch2[subhap.batch2$ID %in% pheno$subjid[pheno$Batch=='II'], ]
subhap <- rbind(subhap.batch1, subhap.batch2)

batch.subhap <- batch.subhap[batch.subhap$TOPMed_ID %in% c(subhap.batch1$ID, subhap.batch2$ID), ]
batch.subhap <- merge(batch.subhap, subhap, by.x = 'TOPMed_ID', by.y = 'ID')

common.ID <- batch.subhap$TOPMed_ID[(batch.subhap$hap1 == batch.subhap$Hap1 & batch.subhap$hap2 == batch.subhap$Hap2) | (batch.subhap$hap1 == batch.subhap$Hap2 & batch.subhap$hap2 == batch.subhap$Hap1)]
common_subset <- pheno$ID %in% common.ID

pred.list   <- 'H1B1G1'
pcs.formula <- sapply(1:10, function(x){paste0('PC',x)})
pcs.formula <- paste(pcs.formula, collapse = ' + ')
euro_subset <- pheno$race %in% 1 & common_subset

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
    effect <- round(result$estimate[result$term == predictor], 3)
    se <- round(result$std.error[result$term == predictor], 3)
    p  <- result$p.value[result$term == predictor]
    freq <- round(sum(pheno[,predictor][assoc_subject])/(2*n), 3)
    output.row <- c('Euro', '1', .pheno, predictor, freq, n, effect, se, p)
    output.row <- unname(output.row)
    model.all.data <- rbind(model.all.data, output.row)
  }
}

# Calc FDR
model.all.data <- model.all.data[2:nrow(model.all.data), ]
model.all.data$Freq >= 0.05
model.all.data$FDR <- NA
model.all.data <- model.all.data[model.all.data$Predictor != 'H2A2G1', ]
# four global phenotypes: CTE, CTEstage, Dementia, AT8 positive tissue
four.pheno   <- c(pheno_list_binary, 'CTEstage', 'AT8PositiveTissue.log')
subhaplotypes <- pred.list[1:9]
dosages       <- pred.list[10:14]
model.all.data$FDR[model.all.data$Phenotype %in% four.pheno & 
                     model.all.data$Predictor %in% subhaplotypes] <- p.adjust(model.all.data$P.value[model.all.data$Phenotype %in% four.pheno & 
                                                                                                       model.all.data$Predictor %in% subhaplotypes],
                                                                              method = 'fdr')

model.all.data$FDR[model.all.data$Phenotype %in% four.pheno & 
                     model.all.data$Predictor %in% dosages] <- p.adjust(model.all.data$P.value[model.all.data$Phenotype %in% four.pheno & 
                                                                                                       model.all.data$Predictor %in% dosages],
                                                                              method = 'fdr')
# eleven endophenotypes
endophenotypes <- pheno_list_categorical[pheno_list_categorical != 'CTEstage']
model.all.data$FDR[model.all.data$Phenotype %in% endophenotypes & 
                     model.all.data$Predictor %in% subhaplotypes] <- p.adjust(model.all.data$P.value[model.all.data$Phenotype %in% endophenotypes & 
                                                                                                       model.all.data$Predictor %in% subhaplotypes],
                                                                              method = 'fdr')
model.all.data$FDR[model.all.data$Phenotype %in% endophenotypes & 
                     model.all.data$Predictor %in% dosages] <- p.adjust(model.all.data$P.value[model.all.data$Phenotype %in% endophenotypes & 
                                                                                                       model.all.data$Predictor %in% dosages],
                                                                              method = 'fdr')
# save result
outputname <- 'Model1.Euro.assoc.txt'
write.table(model.all.data, outputname, col.names = T, row.names = F, sep = '\t', quote = F)
# ==============================================================================================


