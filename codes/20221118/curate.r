setwd('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20221118')
# curate the data
unite <- haven::read_sav('/restricted/projectnb/cte/Challenge_Project_2022/data/pheno/update_11182022/APOEAnalysis.sav')
write.csv(unite, '/restricted/projectnb/cte/Challenge_Project_2022/data/pheno/update_11182022/APOEAnalysis.csv', row.names = F)
# add one subject exlcusively in old phenotype file
add_one <- read.csv('/restricted/projectnb/cte/Challenge_Project_2022/data/pheno/recent/_updates/APOEAnalysis_UNITE_wo_unmatched.csv', as.is = T)
add_one$numberofregions <- NULL

add_one <- add_one[, match(colnames(unite), colnames(add_one))]
unite <- rbind(unite, add_one[nrow(add_one), ])

batch1 <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/data/genetic/1st_batch/_crary_hrc.imputation/CRARY_chr1.hrc.sample',
                     header = F, as.is = T)

unite.batch1 <- unite[unite$pin %in% batch1$V1, ]

batch2 <- read.csv('/restricted/projectnb/cte/Challenge_Project_2022/data/genetic/2nd_batch/ID_study.csv', as.is = T)
unite$subjid[unite$subjid %in% 'K-0487'] <- 'k-0487'
unite$subjid[unite$subjid %in% 'K-0554'] <- 'k-0554'

unite.batch2 <- unite[unite$subjid %in% batch2$Sample_ID, ]

# finding missing subjects in batch2
batch2.missing <- batch2[!batch2$Sample_ID %in% unite.batch2$subjid, ]
batch2.missing.unite <- batch2.missing[batch2.missing$Study == 'UNITE', ]
write.table(batch2.missing.unite, 'missing_data.txt', col.names = T, row.names = F, quote = F, sep = '\t')

unite.merge <- rbind(unite.batch1, unite.batch2)
write.csv(unite.merge, '/restricted/projectnb/cte/Challenge_Project_2022/data/pheno/update_11182022/APOEAnalysis_UNITE_wo_unmatched.csv',
          row.names = F)
