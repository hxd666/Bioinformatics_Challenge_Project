setwd('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20221103')
# Get the sample IDs of UNITE batch I & II

# read IDs for batch I
gwas1 <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/data/genetic/1st_batch/_crary_hrc.imputation/CRARY_hrc_chr17.dose.fam',
                    as.is = T)
gwas1 <- gwas1$V1

# read phenotype file with Cohort information for batch I
pheno1 <- read.csv('/restricted/projectnb/cte/Challenge_Project_2022/data/pheno/old/APOEAnalysiswLEGEND.csv',
                   as.is = T)
# keep UNITE samples only
pheno1 <- pheno1[pheno1$pin %in% gwas1, ]
# keep sample IDs and Race only
pheno1 <- pheno1[, c("pin", "Race_Self_Report")]
colnames(pheno1) <- c('ID', 'Race')
pheno1$Race <- ifelse(pheno1$Race == 'White', 'White', 
                      ifelse(pheno1$Race == 'Black/African American', 'Black', 'Other'))
write.table(pheno1, 'UNITE_batch1.ID.tbx',col.names = T, row.names = F, sep = '\t', quote = F)

# read IDs for batch II
gwas2 <- read.csv('/restricted/projectnb/cte/Challenge_Project_2022/data/genetic/2nd_batch/ID_study.csv',
                  as.is = T)
# keep UNITE only
gwas2 <- gwas2$Sample_ID[gwas2$Study == 'UNITE']
gwas2[gwas2 == 'k-0487'] <- 'K-0487'
gwas2[gwas2 == 'k-0554'] <- 'K-0554'

# read phenotype file with Cohort information for batch II
pheno2 <- read.csv('/restricted/projectnb/cte/Challenge_Project_2022/data/pheno/recent/APOEAnalysis.csv',
                   as.is = T)
ID2 <- data.frame('ID' = gwas2, stringsAsFactors = F)
ID2$Race <- NA
ID2$Race <- pheno2$race[match(ID2$ID, pheno2$subjid)]
ID2$Race <- ifelse(is.na(ID2$Race), 'Unmatched', 
                   ifelse(ID2$Race == 'Black', 'Black', 
                          ifelse(ID2$Race == 'White', 'White',
                                 ifelse(ID2$Race == '', 'Unknown', 'Other'))))
write.table(ID2, 'UNITE_batch2.ID.tbx', col.names = T, row.names = F, sep = '\t', quote = F)