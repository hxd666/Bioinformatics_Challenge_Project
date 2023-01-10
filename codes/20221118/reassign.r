# script to draw PCA plot
batch1.id <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20221107/UNITE_batch1.ID.tbx', header = T, as.is = T)
batch2.id <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20221107/UNITE_batch2.ID.tbx', header = T, as.is = T)
diagnose.id <- read.table('/restricted/projectnb/cte/Challenge_Project_2022/xdhan/scripts/20221107/DIAGNOSE.ID.tbx', header = T, as.is = T)

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



# Function to Calculate the Centroid:
centroid <- function(dataset, race){
  dataset <- dataset[dataset$Race == race, ]
  pc1.mean <- mean(as.numeric(dataset[,'V2']))
  pc2.mean <- mean(as.numeric(dataset[,'V3']))
  return(c(pc1.mean, pc2.mean))
}

distance <- function(dataset, ref.race, study.race){
  center.point <- centroid(dataset, ref.race)
  dataset <- dataset[dataset$Race == study.race, ]
  distance.x <- sapply(dataset$V2, function(x){(as.numeric(x)-center.point[1])^2})
  distance.y <- sapply(dataset$V3, function(y){(as.numeric(y)-center.point[2])^2})
  distance.square <- distance.x + distance.y
  distance.out <- sqrt(distance.square)
  return(distance.out)
}

sd.distance <- function(dataset, race){
  return(sd(distance(dataset, race, race)))
}

check.1outlier <- function(dataset, ref.race, study.race){
  distance.out <- distance(dataset, ref.race, study.race)
  outlier <- distance.out > sd.distance(dataset, ref.race)
  outlier.id <- dataset[dataset$Race == study.race, ]$ID[outlier]
  return(outlier.id)
}

check.2outlier <- function(dataset, ref.race, study.race){
  distance.out <- distance(dataset, ref.race, study.race)
  outlier <- distance.out > 2*sd.distance(dataset, ref.race)
  outlier.id <- dataset[dataset$Race == study.race, ]$ID[outlier]
  return(outlier.id)
}

check.3outlier <- function(dataset, ref.race, study.race){
  distance.out <- distance(dataset, ref.race, study.race)
  outlier <- distance.out > 3*sd.distance(dataset, ref.race)
  outlier.id <- dataset[dataset$Race == study.race, ]$ID[outlier]
  return(outlier.id)
}

check.cluster <- function(dataset, ref.race, study.race){
  distance.out <- distance(dataset, ref.race, study.race)
  outlier <- distance.out < 2*sd.distance(dataset, ref.race)
  outlier.id <- dataset[dataset$Race == study.race, ]$ID[outlier]
  return(outlier.id)
}

# ==============================================================
# For UNITE Batch1
batch1.plot <- rbind(batch1.plot, rep(0, ncol(batch1.plot)))
batch1.plot[nrow(batch1.plot), 22] <- 'Centroid Black'
batch1.plot[nrow(batch1.plot), 2:3] <- centroid(batch1.plot, 'Ref Black')

batch1.plot <- rbind(batch1.plot, rep(0, ncol(batch1.plot)))
batch1.plot[nrow(batch1.plot), 22] <- 'Centroid White'
batch1.plot[nrow(batch1.plot), 2:3] <- centroid(batch1.plot, 'Ref White')
# check Black
ggplot(batch1.plot[batch1.plot$Race %in% c('Ref Black', 'Centroid Black') | batch1.plot$ID %in% check.1outlier(batch1.plot, 'Ref Black', 'Black'), ], aes(x=V2, y =V3)) + geom_point(aes(color=factor(Race)), alpha = 0.8, cex = 2) +
  xlab("PC1") + ylab("PC2") + labs(col="Race", title = 'UNITE Batch1 1SD') +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18))

ggplot(batch1.plot[batch1.plot$Race %in% c('Ref Black', 'Centroid Black') | batch1.plot$ID %in% check.2outlier(batch1.plot, 'Ref Black', 'Black'), ], aes(x=V2, y =V3)) + geom_point(aes(color=factor(Race)), alpha = 0.8, cex = 2) +
  xlab("PC1") + ylab("PC2") + labs(col="Race", title = 'UNITE Batch1 2SD') +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18))

ggplot(batch1.plot[batch1.plot$Race %in% c('Ref Black', 'Centroid Black') | batch1.plot$ID %in% check.3outlier(batch1.plot, 'Ref Black', 'Black'), ], aes(x=V2, y =V3)) + geom_point(aes(color=factor(Race)), alpha = 0.8, cex = 2) +
  xlab("PC1") + ylab("PC2") + labs(col="Race", title = 'UNITE Batch1 3SD') +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18))

# check White
ggplot(batch1.plot[batch1.plot$Race %in% c('Ref White', 'Centroid White') | batch1.plot$ID %in% check.1outlier(batch1.plot, 'Ref White', 'White'), ], aes(x=V2, y =V3)) + geom_point(aes(color=factor(Race)), alpha = 0.8, cex = 2) +
  xlab("PC1") + ylab("PC2") + labs(col="Race", title = 'UNITE Batch1 1SD') +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18))

ggplot(batch1.plot[batch1.plot$Race %in% c('Ref White', 'Centroid White') | batch1.plot$ID %in% check.2outlier(batch1.plot, 'Ref White', 'White'), ], aes(x=V2, y =V3)) + geom_point(aes(color=factor(Race)), alpha = 0.8, cex = 2) +
  xlab("PC1") + ylab("PC2") + labs(col="Race", title = 'UNITE Batch1 2SD') +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18))

ggplot(batch1.plot[batch1.plot$Race %in% c('Ref White', 'Centroid White') | batch1.plot$ID %in% check.3outlier(batch1.plot, 'Ref White', 'White'), ], aes(x=V2, y =V3)) + geom_point(aes(color=factor(Race)), alpha = 0.8, cex = 2) +
  xlab("PC1") + ylab("PC2") + labs(col="Race", title = 'UNITE Batch1 3SD') +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18))
# ==============================================================


# ==============================================================
# For UNITE batch2
batch2.plot <- rbind(batch2.plot, rep(0, ncol(batch2.plot)))
batch2.plot[nrow(batch2.plot), 22] <- 'Centroid Black'
batch2.plot[nrow(batch2.plot), 2:3] <- centroid(batch2.plot, 'Ref Black')

batch2.plot <- rbind(batch2.plot, rep(0, ncol(batch2.plot)))
batch2.plot[nrow(batch2.plot), 22] <- 'Centroid White'
batch2.plot[nrow(batch2.plot), 2:3] <- centroid(batch2.plot, 'Ref White')
# check Black
ggplot(batch2.plot[batch2.plot$Race %in% c('Ref Black', 'Centroid Black') | batch2.plot$ID %in% check.1outlier(batch2.plot, 'Ref Black', 'Black'), ], aes(x=V2, y =V3)) + geom_point(aes(color=factor(Race)), alpha = 0.8, cex = 2) +
  xlab("PC1") + ylab("PC2") + labs(col="Race", title = 'UNITE batch2 1SD') +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18))

ggplot(batch2.plot[batch2.plot$Race %in% c('Ref Black', 'Centroid Black') | batch2.plot$ID %in% check.2outlier(batch2.plot, 'Ref Black', 'Black'), ], aes(x=V2, y =V3)) + geom_point(aes(color=factor(Race)), alpha = 0.8, cex = 2) +
  xlab("PC1") + ylab("PC2") + labs(col="Race", title = 'UNITE batch2 2SD') +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18))

ggplot(batch2.plot[batch2.plot$Race %in% c('Ref Black', 'Centroid Black') | batch2.plot$ID %in% check.3outlier(batch2.plot, 'Ref Black', 'Black'), ], aes(x=V2, y =V3)) + geom_point(aes(color=factor(Race)), alpha = 0.8, cex = 2) +
  xlab("PC1") + ylab("PC2") + labs(col="Race", title = 'UNITE batch2 3SD') +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18))

# check White
ggplot(batch2.plot[batch2.plot$Race %in% c('Ref White', 'Centroid White') | batch2.plot$ID %in% check.1outlier(batch2.plot, 'Ref White', 'White'), ], aes(x=V2, y =V3)) + geom_point(aes(color=factor(Race)), alpha = 0.8, cex = 2) +
  xlab("PC1") + ylab("PC2") + labs(col="Race", title = 'UNITE batch2 1SD') +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18))

ggplot(batch2.plot[batch2.plot$Race %in% c('Ref White', 'Centroid White') | batch2.plot$ID %in% check.2outlier(batch2.plot, 'Ref White', 'White'), ], aes(x=V2, y =V3)) + geom_point(aes(color=factor(Race)), alpha = 0.8, cex = 2) +
  xlab("PC1") + ylab("PC2") + labs(col="Race", title = 'UNITE batch2 2SD') +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18))

ggplot(batch2.plot[batch2.plot$Race %in% c('Ref White', 'Centroid White') | batch2.plot$ID %in% check.3outlier(batch2.plot, 'Ref White', 'White'), ], aes(x=V2, y =V3)) + geom_point(aes(color=factor(Race)), alpha = 0.8, cex = 2) +
  xlab("PC1") + ylab("PC2") + labs(col="Race", title = 'UNITE batch2 3SD') +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18))
# ==============================================================

# ==============================================================
# For DIAGNOSE
diagnose.plot <- rbind(diagnose.plot, rep(0, ncol(diagnose.plot)))
diagnose.plot[nrow(diagnose.plot), 22] <- 'Centroid Black'
diagnose.plot[nrow(diagnose.plot), 2:3] <- centroid(diagnose.plot, 'Ref Black')

diagnose.plot <- rbind(diagnose.plot, rep(0, ncol(diagnose.plot)))
diagnose.plot[nrow(diagnose.plot), 22] <- 'Centroid White'
diagnose.plot[nrow(diagnose.plot), 2:3] <- centroid(diagnose.plot, 'Ref White')
# check Black
ggplot(diagnose.plot[diagnose.plot$Race %in% c('Ref Black', 'Centroid Black') | diagnose.plot$ID %in% check.1outlier(diagnose.plot, 'Ref Black', 'Black'), ], aes(x=V2, y =V3)) + geom_point(aes(color=factor(Race)), alpha = 0.8, cex = 2) +
  xlab("PC1") + ylab("PC2") + labs(col="Race", title = 'DIAGNOSE 1SD') +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18))

ggplot(diagnose.plot[diagnose.plot$Race %in% c('Ref Black', 'Centroid Black') | diagnose.plot$ID %in% check.2outlier(diagnose.plot, 'Ref Black', 'Black'), ], aes(x=V2, y =V3)) + geom_point(aes(color=factor(Race)), alpha = 0.8, cex = 2) +
  xlab("PC1") + ylab("PC2") + labs(col="Race", title = 'DIAGNOSE 2SD') +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18))

ggplot(diagnose.plot[diagnose.plot$Race %in% c('Ref Black', 'Centroid Black') | diagnose.plot$ID %in% check.3outlier(diagnose.plot, 'Ref Black', 'Black'), ], aes(x=V2, y =V3)) + geom_point(aes(color=factor(Race)), alpha = 0.8, cex = 2) +
  xlab("PC1") + ylab("PC2") + labs(col="Race", title = 'DIAGNOSE 3SD') +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18))

# check White
ggplot(diagnose.plot[diagnose.plot$Race %in% c('Ref White', 'Centroid White') | diagnose.plot$ID %in% check.1outlier(diagnose.plot, 'Ref White', 'White'), ], aes(x=V2, y =V3)) + geom_point(aes(color=factor(Race)), alpha = 0.8, cex = 2) +
  xlab("PC1") + ylab("PC2") + labs(col="Race", title = 'DIAGNOSE 1SD') +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18))

ggplot(diagnose.plot[diagnose.plot$Race %in% c('Ref White', 'Centroid White') | diagnose.plot$ID %in% check.2outlier(diagnose.plot, 'Ref White', 'White'), ], aes(x=V2, y =V3)) + geom_point(aes(color=factor(Race)), alpha = 0.8, cex = 2) +
  xlab("PC1") + ylab("PC2") + labs(col="Race", title = 'DIAGNOSE 2SD') +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18))

ggplot(diagnose.plot[diagnose.plot$Race %in% c('Ref White', 'Centroid White') | diagnose.plot$ID %in% check.3outlier(diagnose.plot, 'Ref White', 'White'), ], aes(x=V2, y =V3)) + geom_point(aes(color=factor(Race)), alpha = 0.8, cex = 2) +
  xlab("PC1") + ylab("PC2") + labs(col="Race", title = 'DIAGNOSE 3SD') +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18))
# ==============================================================
