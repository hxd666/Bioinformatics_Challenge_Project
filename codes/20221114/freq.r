setwd('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20221114')

subhap1 <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/imputation/1st_batch/subhaplotypes/1st_batch.topmed.mapt.subhaplotypes.txt',
                      header = T, as.is = T)
subhap2 <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/imputation/2nd_batch/subhaplotypes/2nd_batch.topmed.mapt.subhaplotypes.txt',
                      header = T, as.is = T)

unite1 <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20221107/UNITE_batch1.ID.tbx',
                     header = T, as.is = T)
subhap1$ID <- sapply(subhap1$TOPMed_ID, function(x){strsplit(x,'_')[[1]][1]})
unite1hap <- merge(subhap1, unite1, by = 'ID')

unite2 <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20221107/UNITE_batch2.ID.tbx',
                     header = T, as.is = T)
subhap2$ID <- subhap2$TOPMed_ID
unite2hap <- merge(subhap2, unite2, by = 'ID')


diagnose <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20221107/DIAGNOSE.ID.tbx',
                       header = T, as.is = T)
diagnosehap <- merge(subhap2, diagnose, by = 'ID')


table(c(unite1hap$Hap1, unite1hap$Hap2))/(178*2)

table(c(unite2hap$Hap1, unite2hap$Hap2))/(383*2)

table(c(diagnosehap$Hap1, diagnosehap$Hap2))/(244*2)

table(c(unite1hap$Hap1, unite1hap$Hap2,unite2hap$Hap1, unite2hap$Hap2))/((178+383)*2)

