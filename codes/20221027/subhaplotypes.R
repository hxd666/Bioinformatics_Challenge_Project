# Generate sub-haplotypes based on 12 surrogated markers
rm(list = ls())

study <- commandArgs(T)[1]

haplotypes <- "/restricted/projectnb/cte/Challenge_Project_2022/imputation/XXXX_batch/subhaplotypes/XXXX_batch.topmed.phased.impute2.phased.with.ref.haps.surrogate"
haplotypes <- gsub("XXXX", study, haplotypes)

haplotypes.12 <- read.table(haplotypes, header = F, as.is = T)

haplotypes.sample <- "/restricted/projectnb/cte/Challenge_Project_2022/imputation/XXXX_batch/2nd_phasing/XXXX_batch.topmed.phased.with.ref.sample"
haplotypes.sample <- gsub("XXXX", study, haplotypes.sample)
haplotypes.sample <- read.table(haplotypes.sample, header = T, as.is = T)



haplotypes.subhaplotype <- apply(haplotypes.12[,6:ncol(haplotypes.12)], 2, function(x){paste(x,collapse = "")})

haplotypes.subhaplotype.tmp <- gsub("1","a",haplotypes.subhaplotype)
haplotypes.subhaplotype.tmp <- gsub("0","b",haplotypes.subhaplotype.tmp)

haplotypes.subhaplotype.tmp <- gsub("a","0",haplotypes.subhaplotype.tmp)
haplotypes.subhaplotype.tmp <- gsub("b","1",haplotypes.subhaplotype.tmp)

haplotypes.subhaplotype <- haplotypes.subhaplotype.tmp
rm(haplotypes.subhaplotype.tmp)


haplotypes.subhaplotype.data <- data.frame('TOPMed_ID'=haplotypes.sample$ID_2[2:nrow(haplotypes.sample)],
                                        'Hap1'=haplotypes.subhaplotype[seq(1,length(haplotypes.subhaplotype),2)],
                                        'Hap2'=haplotypes.subhaplotype[seq(2,length(haplotypes.subhaplotype),2)],
                                        stringsAsFactors = F)


convert.table <- data.frame('surrogate'=c('0 1 0 0 0 0 1 0 0 0 0 0',
                                          '0 1 0 0 0 0 1 1 0 0 0 0',
                                          '0 1 0 0 0 0 1 1 1 0 0 0',
                                          '0 1 0 0 0 0 1 1 1 1 0 0',
                                          '0 1 1 0 0 0 1 0 0 0 0 0',
                                          '0 1 1 1 0 0 1 0 0 0 0 0',
                                          '1 0 0 0 1 0 0 0 0 0 1 1',
                                          '1 0 0 0 1 1 0 0 0 0 1 0',
                                          '1 0 0 0 1 1 0 0 0 0 1 1'),
                            'subhaplotype'=c('H1B1G1','H1B1G2','H1B1G3','H1B1G4',
                                             'H1B2G1','H1B3G1','H2A1G2','H2A2G1',
                                             'H2A2G2'))
convert.table$surrogate <- sapply(convert.table$surrogate, function(x){gsub(" ","",x)})

haplotypes.subhaplotype.data$Hap1 <- convert.table$subhaplotype[match(haplotypes.subhaplotype.data$Hap1, convert.table$surrogate)]
haplotypes.subhaplotype.data$Hap2 <- convert.table$subhaplotype[match(haplotypes.subhaplotype.data$Hap2, convert.table$surrogate)]


#haplotypes.subhaplotype.data <- haplotypes.subhaplotype.data[complete.cases(haplotypes.subhaplotype.data), ]


add.subhaplotype <- function(table, subhaplotype){
  for (.single in subhaplotype){
    h1 <- table$Hap1 == .single
    h2 <- table$Hap2 == .single
    sub <- h1 + h2
    table <- cbind(table, sub)
    colnames(table)[ncol(table)] <- .single
  }
  return(table)
}

haplotypes.subhaplotype.data <- add.subhaplotype(haplotypes.subhaplotype.data, 
                                              c('H1B1G1','H1B1G2','H1B1G3','H1B1G4',
                                                'H1B2G1','H1B3G1','H2A1G2','H2A2G1','H2A2G2'))

haplotypes.subhaplotype.data$H1_beta <- rowSums(haplotypes.subhaplotype.data[,c('H1B1G1','H1B1G2','H1B1G3','H1B1G4')])+
  2*haplotypes.subhaplotype.data[,'H1B2G1'] + 3*haplotypes.subhaplotype.data[,'H1B3G1']

haplotypes.subhaplotype.data$H2_alpha <- haplotypes.subhaplotype.data[,'H2A1G2'] + 
  2* rowSums(haplotypes.subhaplotype.data[,c('H2A2G1','H2A2G2')])

haplotypes.subhaplotype.data$H1_gamma <- rowSums(haplotypes.subhaplotype.data[,c('H1B1G1','H1B2G1','H1B3G1')]) +
  2*haplotypes.subhaplotype.data[,'H1B1G2'] + 3*haplotypes.subhaplotype.data[,'H1B1G3'] + 4*haplotypes.subhaplotype.data[,'H1B1G4']

haplotypes.subhaplotype.data$H2_gamma <- 2 * rowSums(haplotypes.subhaplotype.data[,c('H2A1G2','H2A2G2')]) + 
  haplotypes.subhaplotype.data[,'H2A2G1']

haplotypes.dosage.data <- "/restricted/projectnb/cte/Challenge_Project_2022/imputation/XXXX_batch/subhaplotypes/XXXX_batch.topmed.phased.impute2.surrogate.dosage"
haplotypes.dosage.data <- gsub("XXXX", study, haplotypes.dosage.data)
haplotypes.dosage.data <- read.table(haplotypes.dosage.data, header = F, as.is = T)

merge.dosage <- function(x){
  dos <- c()
  for (.i in 1:(length(x)/3)){
    dos <- c(dos, x[.i*3-1] + 2*x[.i*3])
  }
  return(dos)
}

haplotypes.dosage <- apply(haplotypes.dosage.data[,6:ncol(haplotypes.dosage.data)],1,merge.dosage)
haplotypes.dosage.table <- data.frame('TOPMed_ID'=haplotypes.sample$ID_2[2:nrow(haplotypes.sample)],
                                   stringsAsFactors = F)
haplotypes.dosage.table$H1_dosage <- 2 - haplotypes.dosage[,1]
haplotypes.dosage.table$H1_beta_dosage <- rowSums(haplotypes.dosage[,2:4])
haplotypes.dosage.table$H2_alpha_dosage <- rowSums(haplotypes.dosage[,5:6])
haplotypes.dosage.table$H1_gamma_dosage <- rowSums(haplotypes.dosage[,7:10])
haplotypes.dosage.table$H2_gamma_dosage <- rowSums(haplotypes.dosage[,11:12])

haplotypes.final <- merge(haplotypes.subhaplotype.data, haplotypes.dosage.table, by = 'TOPMed_ID')

output.name <- "/restricted/projectnb/cte/Challenge_Project_2022/imputation/XXXX_batch/subhaplotypes/XXXX_batch.topmed.mapt.subhaplotypes.txt"
output.name <- gsub("XXXX", study, output.name)
write.table(haplotypes.final, output.name, col.names = T, row.names = F, quote = F, sep = '\t')



