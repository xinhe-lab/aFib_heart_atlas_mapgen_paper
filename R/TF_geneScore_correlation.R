library(ArchR)
library(ComplexHeatmap)
setwd('/project2/gca/aselewa/heart_atlas_project/')
source('R/analysis_utils.R')
source_overlapCNN()

KNN <- 100

satac <- loadArchRProject('ArchR/ArchR_heart_latest_noAtrium/')
rD <- getReducedDims(ArchRProj = satac, reducedDims = 'harmony', corCutOff = 0.75)
idx <- sample(seq_len(nrow(rD)), 5000, replace = F)

knnObj <- computeKNN(data = rD, query = rD[idx,], k = KNN)
keepKnn <- determineOverlapCpp(knnObj, floor(0.1 * KNN))
knnObj <- knnObj[keepKnn==0,]

meta.cell.id <- rep(0, nCells(satac))
for(i in 1:nrow(knnObj)){
    meta.cell.id[knnObj[i, ]] <- i
}

satac$metaCells <- factor(meta.cell.id)
motif.mat.obj <- getGroupSE(ArchRProj = satac, 
                            useMatrix = "MotifMatrix", 
                            groupBy = "metaCells",
                            divideN = T, 
                            scaleTo = NULL)

motif.mat <- motif.mat.obj@assays@data$MotifMatrix
rowd <- rowData(motif.mat.obj)
motif.mat <- motif.mat[as.character(rowd$seqnames) == "z",]
rownames(motif.mat) <- rowd$name[as.character(rowd$seqnames) == "z"]
motif.mat <- motif.mat[,-1]
rownames(motif.mat) <- sub('_.*','',rownames(motif.mat))
#motif.mat <- sweep(x = motif.mat - rowMeans(motif.mat), MARGIN = 1, STATS = matrixStats::rowSds(motif.mat), FUN = '/')

geneScore.mat.obj <- getGroupSE(ArchRProj = satac, 
                            useMatrix = "GeneScoreMatrix", 
                            groupBy = "metaCells",
                            divideN = T, 
                            scaleTo = 10^5)

geneScore.mat <- log2(geneScore.mat.obj@assays@data$GeneScoreMatrix+1)
rownames(geneScore.mat) <- rowData(geneScore.mat.obj)$name
geneScore.mat <- geneScore.mat[,-1] 
geneScore.mat <- geneScore.mat[rowVars(geneScore.mat) > 0, ]
#geneScore.mat <- sweep(x = geneScore.mat - rowMeans(geneScore.mat), MARGIN = 1, STATS = matrixStats::rowSds(geneScore.mat), FUN = '/')

same.genes <- intersect(rownames(geneScore.mat), rownames(motif.mat))

geneScore.mat <- geneScore.mat[same.genes,]
motif.mat <- motif.mat[same.genes,]

cor.res <- rep(0, length(same.genes))
for(i in 1:length(same.genes)){
    cor.res[i] <- cor(geneScore.mat[i,], motif.mat[i,])
}

gene.score.cor <- data.frame(gene_name = rownames(geneScore.mat), GeneScore_Correlation=cor.res)
gene.score.cor %>% write_tsv('TF_Correlation_Results/All_TFMotif_GeneScore_Correlations_CisBP.tsv')

plot(geneScore.mat["MEF2A",], motif.mat["MEF2A",], xlab='log2 Gene Score', ylab='chromVar Motif Z-score', main='MEF2A')

rna.results <- read_tsv('TF_Correlation_Results/All_TFMotif_Expression_Correlations_CisBP.tsv')
rna.results$geneScore_Cor <- cor.res[match(rna.results$gene_name, same.genes)]

rna.results <- rna.results[!is.na(rna.results$geneScore_Cor),]
rna.results <- rna.results[!duplicated(rna.results$gene_name),]

fit <- lm(rna.results$correlation ~ rna.results$geneScore_Cor)
rna.results$correlation_fitted <- fit$fitted.values

tfs.to.highlight <- c("TBX5","GATA4","MEF2A","TEAD1","HAND2","PRRX1","NKX25","ESRRB","ESRRA","MEF2C","ETV1")
idx <- match(tfs.to.highlight, rna.results$gene_name)
idx <- idx[!is.na(idx)]
rna.results$lab <- ""
rna.results$lab[idx] <- rna.results$gene_name[idx]

pdf('manuscript_figures/figure2/TF_Motif_genescore_vs_RNA.pdf', width=8, height=6)

ggplot(rna.results, aes(x = geneScore_Cor, y = correlation, label = lab)) + 
    geom_point() + 
    ggClean() +
    geom_abline(slope = 1, intercept = 0, color='blue') + 
    xlab('Pearson GeneScore + MetaCells') + 
    ylab('Pearson RNA + Pseudo-bulk')  + 
    geom_line(aes(x = geneScore_Cor, y = correlation_fitted), linetype='dashed', color='red') + ggrepel::geom_text_repel(nudge_x = 0.05, nudge_y=0.05, max.overlaps = 10000)


dev.off()

