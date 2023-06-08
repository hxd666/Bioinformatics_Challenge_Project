# check neural network imputed result
impute.v2.batch1 <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/imputation/1st_batch/subhaplotypes/1st_batch.topmed.mapt.subhaplotypes.txt',
                               header = T, as.is = T)
impute.v2.batch1 <- impute.v2.batch1[, c("TOPMed_ID", "H1B1G1")]

impute.v2.batch2 <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/imputation/2nd_batch/subhaplotypes/2nd_batch.topmed.mapt.subhaplotypes.txt',
                               header = T, as.is = T)
impute.v2.batch2 <- impute.v2.batch2[, c("TOPMed_ID", "H1B1G1")]

impute.v2 <- rbind(impute.v2.batch1, impute.v2.batch2)


impute.NN <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/imputation/neural_net/Prediction_data_updated_PredResult_LATEST.txt',
                        header = T, as.is = T, sep = ',')
impute.NN <- impute.NN[order(impute.NN$ID), ]

impute.NN$hap1.n <- ifelse(impute.NN$hap1 == 'H1B1G1', 1, 0)
impute.NN$hap2.n <- ifelse(impute.NN$hap2 == 'H1B1G1', 1, 0)
impute.NN$H1B1G1 <- impute.NN$hap1.n + impute.NN$hap2.n

nn.id1 <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/_be562/data/pred1.fam',
                     as.is = T)
nn.id2 <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/_be562/data/pred2.fam',
                     as.is = T)
impute.NN$ID <- c(nn.id1$V2, nn.id2$V2)

impute.v2$H1B1G1.nn <- impute.NN$H1B1G1[match(impute.v2$TOPMed_ID, impute.NN$ID)]

impute.v2$H1B1G1[is.na(impute.v2$H1B1G1)] <- 3

sum(impute.v2$H1B1G1 == impute.v2$H1B1G1.nn)


# check without reorder

# check neural network imputed result
impute.v2.batch1 <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/imputation/1st_batch/subhaplotypes/1st_batch.topmed.mapt.subhaplotypes.txt',
                               header = T, as.is = T)
impute.v2.batch1 <- impute.v2.batch1[, c("TOPMed_ID", "H1B1G1")]

impute.v2.batch2 <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/imputation/2nd_batch/subhaplotypes/2nd_batch.topmed.mapt.subhaplotypes.txt',
                               header = T, as.is = T)
impute.v2.batch2 <- impute.v2.batch2[, c("TOPMed_ID", "H1B1G1")]

impute.v2 <- rbind(impute.v2.batch1, impute.v2.batch2)


impute.NN <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/imputation/neural_net/Prediction_data_updated_PredResult_LATEST.txt',
                        header = T, as.is = T, sep = ',')
#impute.NN <- impute.NN[order(impute.NN$ID), ]

impute.NN$hap1.n <- ifelse(impute.NN$hap1 == 'H1B1G1', 1, 0)
impute.NN$hap2.n <- ifelse(impute.NN$hap2 == 'H1B1G1', 1, 0)
impute.NN$H1B1G1 <- impute.NN$hap1.n + impute.NN$hap2.n

nn.id1 <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/_be562/data/pred1.fam',
                     as.is = T)
nn.id2 <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/_be562/data/pred2.fam',
                     as.is = T)
impute.NN$ID <- c(nn.id1$V2, nn.id2$V2)

impute.v2$H1B1G1.nn <- impute.NN$H1B1G1[match(impute.v2$TOPMed_ID, impute.NN$ID)]

impute.v2$H1B1G1[is.na(impute.v2$H1B1G1)] <- 3

sum(impute.v2$H1B1G1 == impute.v2$H1B1G1.nn)
