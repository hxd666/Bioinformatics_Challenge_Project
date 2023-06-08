# Compare the SNPs in two Regions 

alnloci <- read.table('GRCh38_alnloci_mpileup_chrm.snps', header = T, as.is = T)
topmed  <- read.table('TOPMed.CTE.hg38.snps', header = T, as.is = T)

alnloci$LocID <- paste(alnloci$CHR, alnloci$BP, sep = "-")
alnloci <- alnloci[, c("LocID", "ALT", "REF")]
colnames(alnloci) <- c('LocID', 'ALT.aln', 'REF.aln')

topmed$LocID <- paste(topmed$CHR, topmed$BP, sep = "-")
topmed <- topmed[, c("LocID", "ALT", "REF")]
colnames(topmed) <- c('LocID', 'ALT.top', 'REF.top')

ds <- merge(alnloci, topmed, by = 'LocID')

ds_matched <- (toupper(ds$ALT.aln) == toupper(ds$ALT.top) & toupper(ds$REF.aln) == toupper(ds$REF.top)) | (toupper(ds$ALT.aln) == toupper(ds$REF.top) & toupper(ds$REF.aln) == toupper(ds$ALT.top))

ds_diff <- ds[!ds_matched, ]
write.table(ds_diff, 'Umatched.SNPs', row.names = F, col.names = T, quote = F, sep = '\t')
  