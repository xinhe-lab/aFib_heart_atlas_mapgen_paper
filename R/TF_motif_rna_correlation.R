library(ArchR)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
setwd('/project2/gca/aselewa/heart_atlas_project/')
source('R/analysis_utils.R')

ideal_order <- c("Cardiomyocyte","Smooth Muscle","Pericyte","Endothelial","Fibroblast","Neuronal", "Lymphoid","Myeloid")

srna <- readRDS('seurat/Heart_RNA_Processed_Combined_NoAtrium.rds')
satac <- loadArchRProject('ArchR/ArchR_heart_latest_noAtrium/')

motif.mat.obj <- getGroupSE(ArchRProj = satac, 
                            useMatrix = "MotifMatrix", 
                            groupBy = "CellTypes", 
                            divideN = T, 
                            scaleTo = NULL)

normal.exp <- Seurat::AverageExpression(object = srna)

motif.mat <- motif.mat.obj@assays@data$MotifMatrix
rowd <- rowData(motif.mat.obj)
motif.mat <- motif.mat[as.character(rowd$seqnames) == "deviations",]
rownames(motif.mat) <- rowd$name[as.character(rowd$seqnames) == "deviations"]

normal.exp.mat <- log2(as.matrix(normal.exp$RNA) + 1)
normal.exp.mat <- sweep(x = normal.exp.mat - rowMeans(normal.exp.mat), MARGIN = 1, STATS = matrixStats::rowSds(normal.exp.mat), FUN = '/')
motif.enrich.df <- sweep(x = motif.mat - rowMeans(motif.mat), MARGIN = 1, STATS = matrixStats::rowSds(motif.mat), FUN = '/')
rownames(motif.enrich.df) <- sub('_.*', '', rownames(motif.enrich.df))

normal.exp.mat <- normal.exp.mat[,ideal_order]
motif.enrich.df <- motif.enrich.df[,ideal_order]

ct.specific.cor <- list()
motif.mat.list <- list()
for(ct in ideal_order){
    curr.tf <- ct.specific.tf[[ct]]
    same.genes <- intersect(curr.tf, rownames(normal.exp.mat))
    ### CORRELATION
    corr.res <- rep(0, length(same.genes))
    for(i in 1:length(same.genes)){
        corr.res[i] <- cor(as.numeric(normal.exp.mat[same.genes[i],]), as.numeric(motif.enrich.df[same.genes[i],]))
    }
    curr.df <- data.frame(gene_name = same.genes, correlation = corr.res, celltype = ct, stringsAsFactors = F)
    ct.specific.cor[[ct]] <- curr.df
}

ct.specific.cor.df <- Reduce(rbind, ct.specific.cor)
ct.specific.cor.df %>% write_tsv('manuscript_figures/figure2/All_TFMotif_Expression_Correlations_CisBP.tsv')
ct.specific.cor.df.filt <- ct.specific.cor.df[ct.specific.cor.df$correlation > 0.7,]
ct.specific.cor.df.filt <- ct.specific.cor.df.filt[!duplicated(ct.specific.cor.df.filt$gene_name),]

motif.enrich.df.cut <- motif.enrich.df[ct.specific.cor.df.filt$gene_name,]
tf.exp.mat.cut <- normal.exp.mat[ct.specific.cor.df.filt$gene_name,]

nClust = 8
set.seed(100)
row.clust.res <- kmeans(motif.enrich.df.cut, centers = nClust)
labs <- row.clust.res$cluster
rAnno <- rowAnnotation(cluster = factor(sort(labs)))
Heatmap(motif.enrich.df.cut[names(sort(labs)),], cluster_rows = F, cluster_columns = F, left_annotation = rAnno, col = circlize::colorRamp2(c(-2,0,2), c("lightblue","white","firebrick")))

lab_map <- data.frame(num = 1:nClust,
                      celltype = c("Lym1","Neuro", "SmooPeri","SmooPeri2","Mye","Fibro", "Lym2","CM"),
                      stringsAsFactors = F)

labs.named <- plyr::mapvalues(x = labs, from = lab_map$num, to = lab_map$celltype)
labs.named <- factor(labs.named, levels = c("CM","SmooPeri","SmooPeri2","Fibro","Neuro","Lym1","Lym2","Mye"))
labs.order <- order(labs.named)
labs.named <- labs.named[labs.order]

motif.enrich.df.cut <- motif.enrich.df.cut[labs.order, ]
tf.exp.mat.cut <- tf.exp.mat.cut[labs.order, ]
ct.specific.cor.df.filt <- ct.specific.cor.df.filt[labs.order, ]

cAnno <- HeatmapAnnotation(which = "column",
                           cluster_names = anno_mark(at = 1:nrow(ct.specific.cor.df.filt),
                                                     labels = ct.specific.cor.df.filt$gene_name,
                                                     side = "top"),
                           r = anno_points(x = ct.specific.cor.df.filt$correlation, ylim = c(0.6, 1)),
                           show_legend = F, height = unit(3, "in"))

p1 <- Heatmap(matrix = t(motif.enrich.df.cut),
              cluster_rows = F, cluster_columns = F, show_column_names = F, show_row_names = T, 
              col = circlize::colorRamp2(c(-2,0.1,2), c("lightblue","white","firebrick")), 
              show_row_dend = F, show_column_dend = F,
              row_title = 'Motif', 
              #rect_gp = gpar(col = "black", lwd = 0.1),
              column_title = NULL,
              name = "Row-Normalized Motif Access", show_heatmap_legend = F,
              use_raster = T)

p2 <- Heatmap(matrix = t(tf.exp.mat.cut),
              cluster_rows = F, cluster_columns = F, show_column_names = F, show_row_names = T,
              col = circlize::colorRamp2(c(-2,0.1,2), c("lightblue","white","firebrick")), 
              top_annotation = cAnno,  
              # rect_gp = gpar(col = "black", lwd = 0.1),
              row_title = 'Expression', 
              show_heatmap_legend = T,
              column_title = NULL,
              name = "Row-Normalized Expression",
              use_raster = T)

pdf('manuscript_figures/figure2/Fig2D_TFMotif_Expression_Heatmap_bigger.pdf', width=20, height=10)
p2 %v% p1
dev.off()



