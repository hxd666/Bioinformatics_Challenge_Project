# script to draw PCA plot
batch1.id <- read.table('UNITE_batch1.ID.tbx', header = T, as.is = T)
batch2.id <- read.table('UNITE_batch2.ID.tbx', header = T, as.is = T)

batch1.pca <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/1st_batch/PCA/unite1.pca.txt',
                         as.is = T)
batch2.pca <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE/2nd_batch/PCA/unite2.pca.txt',
                         as.is = T)

colnames(batch1.pca)[1] <- 'ID'; colnames(batch2.pca)[1] <- 'ID'

batch1.plot <- merge(batch1.id, batch1.pca, by = 'ID')
batch2.plot <- merge(batch2.id, batch2.pca, by = 'ID')

library(ggplot2)

ggplot(batch1.plot, aes(x=V2, y =V3)) + geom_point(aes(color=factor(Race)), alpha = 0.8, cex = 2) +
  xlab("PC1") + ylab("PC2") + labs(col="Race", title = 'UNITE Batch1') +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18))

ggplot(batch1.plot, aes(x=V2, y =V4)) + geom_point(aes(color=factor(Race)), alpha = 0.8, cex = 2) +
  xlab("PC1") + ylab("PC3") + labs(col="Race", title = 'UNITE Batch1') +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18))

ggplot(batch1.plot, aes(x=V3, y =V4)) + geom_point(aes(color=factor(Race)), alpha = 0.8, cex = 2) +
  xlab("PC2") + ylab("PC3") + labs(col="Race", title = 'UNITE Batch1') +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18))

ggplot(batch2.plot, aes(x=V2, y =V3)) + geom_point(aes(color=factor(Race)), alpha = 0.8, cex = 2) +
  xlab("PC1") + ylab("PC2") + labs(col="Race", title = 'UNITE Batch2') +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18))

ggplot(batch2.plot, aes(x=V2, y =V4)) + geom_point(aes(color=factor(Race)), alpha = 0.8, cex = 2) +
  xlab("PC1") + ylab("PC3") + labs(col="Race", title = 'UNITE Batch2') +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18))

ggplot(batch2.plot, aes(x=V3, y =V4)) + geom_point(aes(color=factor(Race)), alpha = 0.8, cex = 2) +
  xlab("PC2") + ylab("PC3") + labs(col="Race", title = 'UNITE Batch2') +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18))
