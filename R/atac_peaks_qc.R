library(ArchR)
library(tidyverse)

setwd('/project2/gca/aselewa/heart_atlas_project/')
source('R/analysis_utils.R')
palet <- readRDS('notebooks/palette.rds')
satac <- loadArchRProject(path = 'ArchR/ArchR_heart_latest_noAtrium/')

mean.peak <- getGroupSE(ArchRProj = satac, 
                        useMatrix = 'PeakMatrix', 
                        groupBy = 'CellTypes', 
                        divideN = T, 
                        scaleTo = 10^5)

mean.peak.mat <- as.data.frame(mean.peak@assays@data$PeakMatrix)
quants <- seq(0.01, 1, 0.05)
mean.peak.quants <- as.data.frame(matrix(0, nrow = length(quants), ncol = ncol(mean.peak.mat)))

celltypes <- colnames(mean.peak.mat)
colnames(mean.peak.quants) <- celltypes
for(ct in celltypes){
    mean.peak.quants[,ct] <- quantile(mean.peak.mat[,ct], quants)
}

mean.peak.quants$quants <- quants
mean.peak.quants <- pivot_longer(data = mean.peak.quants, cols = celltypes)

pdf('manuscript_figures/mean_access_vs_quant.pdf', width=8, height=6)
ggplot(mean.peak.quants, aes(x = quants, y = value, color=name)) + geom_line(size=1) + ggClean() + xlab('Quantiles') + ylab('Mean Accessibility') +
    scale_color_manual(values = palet)
dev.off()


satac$agg <- "agg"
mean.peak.agg <- getGroupSE(ArchRProj = satac, 
                        useMatrix = 'PeakMatrix', 
                        groupBy = 'agg', 
                        divideN = T, 
                        scaleTo = 10^5)

peak.mean <- data.frame(mean.peak.agg@assays@data$PeakMatrix)
peak.mean$peakType <- satac@peakSet$peakType
peaktype.med <- peak.mean %>% group_by(peakType) %>% summarise(medianAcc = mean(agg))

peak.mean.cm <- as.data.frame(mean.peak.mat[,"Cardiomyocyte",drop=FALSE])
peak.mean.cm$peakType <- satac@peakSet$peakType
peaktype.med.cm <- peak.mean.cm %>% group_by(peakType) %>% summarise(medianAcc = mean(Cardiomyocyte))

peaktype.med <- bind_rows(peaktype.med, peaktype.med.cm)
peaktype.med$cellType <- c(rep("Aggregate",4), rep("CM",4))

ggplot(peaktype.med, aes(x = peakType, y = medianAcc, fill=cellType)) + 
    geom_bar(stat='identity', position='dodge') + 
    ggClean() + 
    scale_fill_brewer(palette = 'Set2') + 
    xlab('') +
    ylab('Median')

