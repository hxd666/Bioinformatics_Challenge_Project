# sort data
i <- commandArgs(T)[1]
#i <- '1'

output.fn <- '/restricted/projectnb/cte/Challenge_Project_2022/xdhan/_assoc/20221118/modelXXX.euro.txt'
output.fn <- gsub('XXX', i, output.fn)
output <- read.table(output.fn, header = T, as.is = T)

output1 <- output[output$Batch == 'I', ]
output2 <- output[output$Batch == 'II', ]

metal.fn <- '/restricted/projectnb/cte/Challenge_Project_2022/xdhan/_assoc/20221121/output/modelXXX.euro.metal.1.txt'
metal.fn <- gsub('XXX', i, metal.fn)
metal <- read.table(metal.fn, header = T, as.is = T, sep = '\t')
metal <- metal[match(output1$Phenotype, metal$MarkerName), ]

output1.stats <- output1[, 4:7]
output2.stats <- output2[, 4:7] 
metal.stats   <- metal[, c(12, 4:6)]

output.all <- cbind(output1.stats, output2.stats)
output.all <- cbind(output.all, metal.stats)
output.all <- round(output.all, digits = 4)
output.fn  <- 'modelXXX.stats.txt'
output.fn  <- gsub('XXX', i, output.fn)
write.table(output.all, output.fn, col.names = T, row.names = F, quote = F, sep = '\t')
