# Keep the UNITE GWAS samples without unmatched samples
pheno.new <- read.csv('/restricted/projectnb/cte/Challenge_Project_2022/data/pheno/recent/APOEAnalysis.csv',
                      as.is = T)
pheno.old <- read.csv('/restricted/projectnb/cte/Challenge_Project_2022/data/pheno/old/APOEAnalysiswLEGEND.csv',
                      as.is = T)
pheno.old <- pheno.old[pheno.old$subjid == 'SLI-48', ]
# match the old columns to new columns
pheno.old.colmatch <- pheno.old[, colnames(pheno.old) %in% colnames(pheno.new)]
# unmatched columns
unmatch.col <- colnames(pheno.new)[!(colnames(pheno.new) %in% colnames(pheno.old.colmatch))]
print(unmatch.col)
unmatch.col.old <- colnames(pheno.old)[!(colnames(pheno.old) %in% colnames(pheno.new))]
# add columns one by one
table(pheno.new$race)
pheno.old.colmatch$race <- pheno.old$Race_Self_Report
table(pheno.new$agecogsx)
pheno.old.colmatch$agecogsx <- NA
pheno.old.colmatch$cogsx <- NA
pheno.old.colmatch$fhyn <- NA
pheno.old.colmatch$ageplaypro <- pheno.old$YrsFootballSemiPro
pheno.old.colmatch$numregions <- pheno.old$numberofregions
pheno.old.colmatch$AT8sulcus <- pheno.old$sulcratio
pheno.old.colmatch$AT8crest  <- pheno.old$crestratio
pheno.old.colmatch$locratio <- pheno.old$locusratio
pheno.old.colmatch$CA1ratio <- pheno.old$ca1ratio
pheno.old.colmatch$CA23ratio <- pheno.old$ca23ratio
pheno.old.colmatch$CA4ratio  <- pheno.old$ca4ratio
pheno.old.colmatch$pinid     <- pheno.old$pin
pheno.old.colmatch$npftdtdp  <- NA
pheno.old.colmatch$FTLDtau <- NA
pheno.old.colmatch$midfront <- pheno.old$MidFront
pheno.old.colmatch$inforbfront <- pheno.old$InfOrbFront
pheno.old.colmatch$suptemp <- pheno.old$SupTemp
pheno.old.colmatch$infpar <- pheno.old$InfPar
pheno.old.colmatch$entorhinal <- pheno.old$Entorhinal
pheno.old.colmatch$amygdala <- pheno.old$Amygdala
pheno.old.colmatch$subnigra <- pheno.old$SubNigra
pheno.old.colmatch$loccoer  <- pheno.old$LocCoer
pheno.old.colmatch$agedeath <- pheno.old$AgeDeath

# add SLI-48 to new phenotype file
pheno.merge <- rbind(pheno.new, pheno.old.colmatch)

# gwas for two batches
batch1 <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/data/genetic/1st_batch/_crary_hrc.imputation/CRARY_hrc_chr17.dose.fam', as.is = T)
batch2 <- read.csv('/restricted/projectnb/cte/Challenge_Project_2022/data/genetic/2nd_batch/ID_study.csv', as.is = T)
batch2$Sample_ID[batch2$Sample_ID == 'k-0487'] <- 'K-0487'
batch2$Sample_ID[batch2$Sample_ID == 'k-0554'] <- 'K-0554'
pheno.merge <- pheno.merge[pheno.merge$pin %in% batch1$V1 | pheno.merge$subjid %in% batch2$Sample_ID, ]

write.csv(pheno.merge, '/restricted/projectnb/cte/Challenge_Project_2022/data/pheno/recent/_updates/APOEAnalysis_UNITE_wo_unmatched.csv',
          row.names = F)
