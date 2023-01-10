# script to draw PCA plot
batch1.id <- read.table('unite1.race.txt', header = T, as.is = T)
batch2.id <- read.table('unite2.race.txt', header = T, as.is = T)
diagnose.id <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20221107/DIAGNOSE.ID.tbx', header = T, as.is = T)
diagnose.id$Race[diagnose.id$Race == 'Unmatched'] <- 'Missing'

batch1.pca <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/UNITE_batch1/PCA/UNITE_batch1.pca.txt',
                         as.is = T)
batch2.pca <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/UNITE_batch2/PCA/UNITE_batch2.pca.txt',
                         as.is = T)
diagnose.pca <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/DIAGNOSE/PCA/DIAGNOSE.pca.txt',
                           as.is = T)

euro.id <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/reference/input2.pedind',
                      as.is = T)
afrc.id <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/_PCA/UNITE_DIAGNOSE/reference/input3.pedind',
                      as.is = T)
colnames(euro.id)[1] <- 'ID'; colnames(afrc.id)[1] <- 'ID'

colnames(batch1.pca)[1] <- 'ID'; colnames(batch2.pca)[1] <- 'ID'; colnames(diagnose.pca)[1] <- 'ID'

batch1.plot <- batch1.pca
batch1.plot$Race <- NA
batch1.plot$Race[batch1.plot$ID %in% batch1.id$ID] <- batch1.id$Race[match(batch1.plot$ID[batch1.plot$ID %in% batch1.id$ID], batch1.id$ID)]
batch1.plot$Race[batch1.plot$ID %in% euro.id$ID] <- 'Ref White'
batch1.plot$Race[batch1.plot$ID %in% afrc.id$ID] <- 'Ref Black'

batch2.plot <- batch2.pca
batch2.plot$Race <- NA
batch2.plot$Race[batch2.plot$ID %in% batch2.id$ID] <- batch2.id$Race[match(batch2.plot$ID[batch2.plot$ID %in% batch2.id$ID], batch2.id$ID)]
batch2.plot$Race[batch2.plot$ID %in% euro.id$ID] <- 'Ref White'
batch2.plot$Race[batch2.plot$ID %in% afrc.id$ID] <- 'Ref Black'

diagnose.plot <- diagnose.pca
diagnose.plot$Race <- NA
diagnose.plot$Race[diagnose.plot$ID %in% diagnose.id$ID] <- diagnose.id$Race[match(diagnose.plot$ID[diagnose.plot$ID %in% diagnose.id$ID], diagnose.id$ID)]
diagnose.plot$Race[diagnose.plot$ID %in% euro.id$ID] <- 'Ref White'
diagnose.plot$Race[diagnose.plot$ID %in% afrc.id$ID] <- 'Ref Black'

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

ggplot(diagnose.plot, aes(x=V2, y =V3)) + geom_point(aes(color=factor(Race)), alpha = 0.8, cex = 2) +
  xlab("PC1") + ylab("PC2") + labs(col="Race", title = 'DIAGNOSE') +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18))

ggplot(diagnose.plot, aes(x=V2, y =V4)) + geom_point(aes(color=factor(Race)), alpha = 0.8, cex = 2) +
  xlab("PC1") + ylab("PC3") + labs(col="Race", title = 'DIAGNOSE') +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18))

ggplot(diagnose.plot, aes(x=V3, y =V4)) + geom_point(aes(color=factor(Race)), alpha = 0.8, cex = 2) +
  xlab("PC2") + ylab("PC3") + labs(col="Race", title = 'UNITE diagnose') +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18))


# other races only
ggplot(batch1.plot[batch1.plot$Race %in% c('Other','Unknown','Unmatched'),], aes(x=V2, y =V3)) + geom_point(aes(color=factor(Race)), alpha = 0.8, cex = 2) +
  xlab("PC1") + ylab("PC2") + labs(col="Race", title = 'UNITE Batch1') +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18))
ggplot(batch1.plot[batch1.plot$Race=='Other',], aes(x=V3, y =V4)) + geom_point(aes(color=factor(Race)), alpha = 0.8, cex = 2) +
  xlab("PC1") + ylab("PC3") + labs(col="Race", title = 'UNITE Batch1') +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18))

ggplot(batch2.plot[batch2.plot$Race %in% c('Other','Unknown'),], aes(x=V2, y =V3)) + geom_point(aes(color=factor(Race)), alpha = 0.8, cex = 1) +
  xlab("PC1") + ylab("PC2") + labs(col="Race", title = 'UNITE Batch2') +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18))

ggplot(diagnose.plot[diagnose.plot$Race %in% c('Other','Unknown'),], aes(x=V2, y =V3)) + geom_point(aes(color=factor(Race)), alpha = 0.8, cex = 1.5) +
  xlab("PC1") + ylab("PC2") + labs(col="Race", title = 'DIAGNOSE') +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18))


