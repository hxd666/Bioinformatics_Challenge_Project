# conver the best guesses to 1 and non-best guesses to 0

rm(list = ls())
study <- commandArgs(T)[1]

surrogate <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/data/genetic/reference/mapt/surrogate.12.pos.txt', header = F, as.is = T)

impute2 <- "/restricted/projectnb/cte/Challenge_Project_2022/imputation/XXXX_batch/impute/XXXX_batch.topmed.phased.impute2"
impute2 <- gsub("XXXX", study, impute2)

impute2 <- read.table(impute2, header = F, as.is = T)

impute2.surrogate <- impute2[impute2$V3 %in% surrogate[,3], ]

biggestone <- function(x){
  max_ind <- c()
  for (.i in 1:(length(x)/3)){
    .max <- which.max(c(x[ 3*.i-2 ], x[ 3*.i-1 ], x[ 3*.i ]))
    .max <- .max + (.i-1)*3
    max_ind <- c(max_ind, .max)
  }
  return(max_ind)
}

convert_guesses <- function(x){
  x.ind <- biggestone(x)
  x[x.ind] <- 1
  x[setdiff(1:length(x), x.ind)] <- 0
  return(x)
}

impute2.surrogate[,6:ncol(impute2.surrogate)] <- t(apply(impute2.surrogate[,6:ncol(impute2.surrogate)],1,convert_guesses))

impute2[impute2$V3 %in% surrogate[,3],] <- impute2.surrogate

impute2$V1 <- "17"

impute2$V4[impute2$V3 %in% surrogate$V3] <- "1"
impute2$V5[impute2$V3 %in% surrogate$V3] <- "2"

outfile <- "/restricted/projectnb/cte/Challenge_Project_2022/imputation/XXXX_batch/2nd_phasing/XXXX_batch.topmed.phased.impute2"
outfile <- gsub("XXXX", study, outfile)
write.table(impute2, outfile, row.names = F, col.names = F, 
            quote = F, sep = ' ')



