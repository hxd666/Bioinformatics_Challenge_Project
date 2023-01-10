setwd('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20221103')
# Get the sample IDs of DIAGNOSE

# read IDs for batch II
gwas2 <- read.csv('/restricted/projectnb/cte/Challenge_Project_2022/data/genetic/2nd_batch/ID_study.csv',
                  as.is = T)
# keep UNITE only
gwas2 <- gwas2$Sample_ID[gwas2$Study == 'DIAGNOSE']


# read phenotype file with Cohort information for batch II
pheno2 <- read.csv('/restricted/projectnb/cte/Challenge_Project_2022/data/pheno/recent/custom_merged_20220331.csv',
                   as.is = T)
ID2 <- data.frame('ID' = gwas2, stringsAsFactors = F)
ID2$Race <- NA
ID2$Race <- pheno2$racecat_combined[match(ID2$ID, pheno2$subject_id)]
ID2$Race <- ifelse(is.na(ID2$Race), 'Unmatched', 
                   ifelse(ID2$Race == 3, 'Black', 
                          ifelse(ID2$Race == 5, 'White',
                                 ifelse(ID2$Race == 0, 'Unknown', 'Other'))))
write.table(ID2, 'DIAGNOSE.ID.tbx', col.names = T, row.names = F, sep = '\t', quote = F)
