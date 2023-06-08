# sort the metal results
suppressPackageStartupMessages(library(dplyr))
metal <- read.table('UNITE_H1B1G1.assoc.1.txt', header = T, as.is = T)
metal <- metal %>% filter(HetDf >= 1)

template <- read.table('../results/Model1.Euro.assoc.batch2.txt', header = T, as.is = T)
template <- template %>% filter(SNP %in% metal$MarkerName)

metal <- metal[match(template$SNP, metal$MarkerName), ]
all(metal$MarkerName == template$SNP)

template$Freq <- metal$Freq1
template$Effect <- metal$Effect
template$N    <- metal$N
template$SE   <- metal$StdErr
template$P.value <- metal$P.value

template$A1 <- NULL; template$A2 <- NULL; template$SNP <- NULL; template$FDR <- NULL

load('../imputatoin.RData')
template$FDR
four.pheno   <- c(pheno_list_binary, 'CTEstage', 'AT8PositiveTissue.log')
pred.list <- read.table('../predictors_list', header = F, as.is = T)
pred.list <- pred.list$V1
subhaplotypes <- pred.list[1:9]
dosages       <- pred.list[10:14]
template$FDR[template$Phenotype %in% four.pheno & 
                     template$Predictor %in% subhaplotypes] <- p.adjust(template$P.value[template$Phenotype %in% four.pheno & 
                                                                                                       template$Predictor %in% subhaplotypes],
                                                                              method = 'fdr')

template$FDR[template$Phenotype %in% four.pheno & 
                     template$Predictor %in% dosages] <- p.adjust(template$P.value[template$Phenotype %in% four.pheno & 
                                                                                                 template$Predictor %in% dosages],
                                                                        method = 'fdr')
# eleven endophenotypes
endophenotypes <- pheno_list_categorical[pheno_list_categorical != 'CTEstage']
template$FDR[template$Phenotype %in% endophenotypes & 
                     template$Predictor %in% subhaplotypes] <- p.adjust(template$P.value[template$Phenotype %in% endophenotypes & 
                                                                                                       template$Predictor %in% subhaplotypes],
                                                                              method = 'fdr')
template$FDR[template$Phenotype %in% endophenotypes & 
                     template$Predictor %in% dosages] <- p.adjust(template$P.value[template$Phenotype %in% endophenotypes & 
                                                                                                 template$Predictor %in% dosages],
                                                                        method = 'fdr')

write.table(template, '../results/Model.Euro.assoc.meta.txt', col.names = T, row.names = F, sep = '\t',
            quote = F)
