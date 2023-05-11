library(ArchR)
library(tidyverse)
library(ComplexHeatmap)
setwd('/project2/gca/aselewa/heart_atlas_project/')
source('R/analysis_utils.R')
palette <- readRDS('notebooks/palette.rds')

satac <- suppressMessages(loadArchRProject('ArchR/ArchR_heart_latest_noAtrium/'))
peak.markers <- readRDS('ArchR/ArchR_heart_latest_noAtrium/PeakCalls/DA_MARKERS_FDRP_1_log2FC_1.rds')
celltype_ideal_order <- c("Cardiomyocyte","Smooth Muscle","Pericyte","Endothelial","Fibroblast","Neuronal", "Lymphoid","Myeloid")
peak.set <- satac@peakSet
peak.set$peakID <- GRToString(peak.set)

srna <- readRDS('seurat/Heart_RNA_Processed_Combined_NoAtrium.rds')
srna.exp <- Seurat::AverageExpression(srna)
srna.exp.mat <- as.matrix(log2(srna.exp$RNA + 1))
srna.exp.mat.z <- sweep(srna.exp.mat - rowMeans(srna.exp.mat), MARGIN = 1, STATS = matrixStats::rowSds(srna.exp.mat), FUN = '/')

finemap.res <- readr::read_tsv('eQTL_enrich/gtex_finemapping/GTEx_v8_finemapping_DAPG/Heart_LV_Finemapping_CS95.txt', col_names = F)
colnames(finemap.res) <- c("tissue","gene","cluster_id","cluster_pip","variant_id","variant_pip")
finemap.res <- finemap.res[finemap.res$variant_id %in% eqtl.snp.gr$SNP,] %>% 
    dplyr::select(variant_id, gene) %>% rename(gene_id = gene, SNP = variant_id) %>% mutate(gene_id = sub('\\..*', '', gene_id))

ensembl.to.symbol <- AnnotationDbi::select(x = org.Hs.eg.db::org.Hs.eg.db, keys = finemap.res$gene_id, columns = "SYMBOL", keytype = "ENSEMBL") %>% 
    rename(gene_id = ENSEMBL) %>%
    group_by(gene_id) %>% slice(1)

finemap.res <- left_join(finemap.res, ensembl.to.symbol)

peak.markers.union <- peak.markers
eqtl.overlap.peaks <- lapply(peak.markers.union, function(x){plyranges::join_overlap_inner(eqtl.snp.gr, x)})
eqtl.overlap.peaks$union <- NULL
ns <- rep(names(eqtl.overlap.peaks), lengths(eqtl.overlap.peaks))
eqtl.overlap.peaks <- unlist(GRangesList(eqtl.overlap.peaks))
eqtl.overlap.peaks$cellType <- factor(ns, levels = celltype_ideal_order)
eqtl.overlap.peaks.df <- eqtl.overlap.peaks %>% as_tibble() %>% select(SNP, cellType)

finemap.res <- finemap.res %>% left_join(., eqtl.overlap.peaks.df, on='SNP')

finemap.res.inRNA <- finemap.res[finemap.res$SYMBOL %in% rownames(srna.exp.mat),] %>% arrange(cellType) %>% filter(!is.na(cellType))
srna.egene.mat <- srna.exp.mat.z[finemap.res.inRNA$SYMBOL, ]

mat.list <- list()
for(ct in celltype_ideal_order){
    curr.mat <- srna.egene.mat[finemap.res.inRNA$cellType == ct,]
    mat.list[[ct]] <- curr.mat[order(curr.mat[, ct], decreasing = T),]
}
srna.egene.mat.sorted <- Reduce(rbind, mat.list)
srna.egene.mat.sorted <- srna.egene.mat.sorted[,celltype_ideal_order]
rownames(srna.egene.mat.sorted) <- make.unique(rownames(srna.egene.mat.sorted))

genes.to.show <- c("ESRRB","MAST4","SOX9","FGD4","CAMK1D","COL18A1","EVA1B","DLK1","IL32.1","CDK11B","MICB.1","STX3")

pdf('manuscript_figures/figure4/Fig4_celltype_specific_eGene_exp.pdf', width=10, height=6)
p <- Heatmap(t(srna.egene.mat.sorted), show_row_names = T, show_column_names = F,
             cluster_rows = F, cluster_columns = F, 
             bottom_annotation  = HeatmapAnnotation(which = "column", 
                                                    cluster = finemap.res.inRNA$cellType,
                                                    col = list(cluster = palette[celltype_ideal_order])),
             top_annotation = HeatmapAnnotation(which = "column", 
                                                genestoshow = anno_mark(at = match(genes.to.show, rownames(srna.egene.mat.sorted)), 
                                                                        labels = sub('\\..*','',genes.to.show), side = 'top')),
             col = circlize::colorRamp2(c(-3, 0, 3), c("lightblue","white","firebrick")),
             heatmap_legend_param = list(direction = "horizontal"))
draw(p, heatmap_legend_side = "bottom")
dev.off()