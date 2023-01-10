setwd('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20221118')

unite.new <- read.csv('/restricted/projectnb/cte/Challenge_Project_2022/data/pheno/update_11182022/APOEAnalysis_UNITE_wo_unmatched.csv',
                      as.is = T)
batch1 <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/UNITE_batch1/unite1.fam',
                     as.is = T)
batch1$Race <- unite.new$race[match(batch1$V1, unite.new$pin)]
batch1 <- batch1[, c("V1", "Race")]
colnames(batch1)[1] <- 'ID'
table(batch1$Race, useNA = 'always')
batch1$Race[batch1$Race == 'White'] <- '1'
batch1$Race <- ifelse(batch1$Race == '1', 'White', 
                      ifelse(batch1$Race == '2', 'Black', 'Other'))
table(batch1$Race, useNA = 'always')
write.table(batch1, 'unite1.race.txt', col.names = T, row.names = F, sep = '\t', quote = F)

batch2 <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/UNITE_batch2/unite2.fam',
                     as.is = T)
batch2$Race <- NA
batch2 <- batch2[, c("V2", "Race")]
colnames(batch2)[1] <- 'ID'
batch2$Race <- unite.new$race[match(batch2$ID, unite.new$subjid)]
batch2$Race[!batch2$ID %in% unite.new$subjid] <- 'Missing'
table(batch2$Race, useNA = 'always')

batch2$Race <- ifelse(batch2$Race == '1', 'White', 
                      ifelse(batch2$Race == '2', 'Black',
                             ifelse(batch2$Race %in% c('4','5','6'), 'Other',
                                    ifelse(is.na(batch2$Race), 'Unknown', batch2$Race))))
batch2$Race[is.na(batch2$Race)] <- 'Unknown'
table(batch2$Race, useNA = 'always')
write.table(batch2, 'unite2.race.txt', col.names = T, row.names = F, sep = '\t', quote = F)
