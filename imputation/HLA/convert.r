# convert dosage data to the format consistent with impute v2 results
subhap <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/imputation/2nd_batch/subhaplotypes/2nd_batch.topmed.mapt.subhaplotypes.txt',
                     header = T, as.is = T)

hla <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/imputation/HLA/2nd_batch/HLA.dosage.txt', header = F, as.is = T)
hla.list <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/imputation/HLA/2nd_batch/2nd_batch.list', header = F, as.is = T)
hla.list <- hla.list$V1

hla.types <- hla$V3
hla <- t(hla[, 6:793])
hla <- as.data.frame(hla, stringsAsFactors = F)
colnames(hla) <- hla.types
hla$TOPMed_ID <- hla.list

subhap <- merge(subhap, hla)
write.table(subhap, '/restricted/projectnb/cte/Challenge_Project_2022/imputation/HLA/2nd_batch/2nd_batch.hla.types.txt',
            col.names = T,
            row.names = F,
            sep = '\t',
            quote = F)
