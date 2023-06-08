# Clean data for RNA-seq analysis
meta <- read.csv('/restricted/projectnb/cte/Challenge_Project_2022/data/RNA/all_cte_meta_nmf.csv',
                 as.is = T)

# Modify IDs in SampleName (change K-XXX to K-0XXX):
samplename <- meta$SampleName
is.K_ <- grepl('K-', samplename) & !(grepl('K-0', samplename))
samplename[is.K_] <- sapply(samplename[is.K_], function(x){gsub('K-','K-0',x)})

meta$newSampleName <- samplename

# load phenotype file
load('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20230201/imputatoin.RData')
sum(pheno$subjid %in% meta$newSampleName)

pheno_sub <- pheno[, c(3, 23:31)]

meta <- merge(meta, pheno_sub, by.x = 'newSampleName', by.y = 'subjid', all.x = T)

rs_1st <- read.table('1st_batch.rs6910507.GT', as.is = T)
rs_1st$V1 <- NULL; rs_1st$V2 <- NULL
rs_1st_list <- read.table('1st_sample', as.is = T)
rs_1st <- cbind(rs_1st_list, t(rs_1st))
colnames(rs_1st) <- c('ID', 'rs6910507')
rs_1st$ID <- sapply(rs_1st$ID, function(x){strsplit(x,'_')[[1]][1]})


rs_2nd <- read.table('2nd_batch.rs6910507.GT', as.is = T)
rs_2nd$V1 <- NULL; rs_2nd$V2 <- NULL
rs_2nd_list <- read.table('2nd_sample', as.is = T)
rs_2nd <- cbind(rs_2nd_list, t(rs_2nd))
colnames(rs_2nd) <- c('ID', 'rs6910507')

rs <- rbind(rs_1st, rs_2nd)
rs$rs6910507 <- ifelse(rs$rs6910507 %in% c('1|0', '1|1', '0|1'), 1, 0)
rs$newSampleName <- pheno$subjid[match(rs$ID, pheno$pin)]
rs$newSampleName[is.na(rs$newSampleName)] <- rs$ID[is.na(rs$newSampleName)]
rs$newSampleName[rs$newSampleName == 'k-0487'] <- 'K-0487'
rs$newSampleName[rs$newSampleName == 'k-0554'] <- 'K-0554'

meta$rs6910507 <- rs$rs6910507[match(meta$newSampleName, rs$newSampleName)]

output.filename <- '/restricted/projectnb/cte/Challenge_Project_2022/data/RNA/all_cte_meta_nmf_cleaned.csv'
write.csv(meta, output.filename, row.names = F, quote = F, col.names = T)
