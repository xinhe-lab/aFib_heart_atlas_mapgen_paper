source('/project2/gca/aselewa/heart_atlas_project/R/analysis_utils.R')
library(Seurat)
library(ArchR)

## Run eqtl_enrichment.Rmd notebook first to generate the files required here!

# load projects
satac<- loadArchRProject('ArchR/ArchR_heart_latest_noAtrium/')
srna <- readRDS('seurat/Heart_RNA_Processed_Combined_NoAtrium.rds')

# average RNA matrix per celltype
srna.exp <- AverageExpression(srna)
srna.exp.mat <- as.matrix(log2(srna.exp$RNA + 1))

# average RNA bulk matrix
Idents(srna) <- "agg"
agg.exp <- AverageExpression(srna)
agg.exp.mat <- as.matrix(log2(agg.exp$RNA + 1))

# partitioning of peaks into DA and shared
peaks.list <- readRDS('ArchR/ArchR_heart_latest_noAtrium/PeakCalls/DA_Plus_Quantiles_peaks_list.rds')

# load eQTL data
eqtl.snp.gr <- readRDS('eQTL_enrich/gtex_finemapping/HighPIP_eQTLs_w_LFSR.gr.rds')
eqtl.snp.gr$gene <- sub('[.].*','',eqtl.snp.gr$gene)

# convert ensembl IDs to symbol

ensembl.to.symbol <- function(ids) {
  AnnotationDbi::select(x = org.Hs.eg.db::org.Hs.eg.db, keys = ids, columns = "SYMBOL", keytype = "ENSEMBL") %>% 
  rename(gene_id = ENSEMBL) %>%
  group_by(gene_id) %>% slice(1) %>%
  filter(!is.na(SYMBOL)) %>% .$SYMBOL
}

# get CM specific eQTLs and shared eQTLs
cm.eqtl.snp.gr <- subsetByOverlaps(eqtl.snp.gr,peaks.list$Cardiomyocyte)
shared.eqtl.snp.gr <- subsetByOverlaps(eqtl.snp.gr,peaks.list$`CM_4+_Shared`)

cm.specific.egenes <- intersect(ensembl.to.symbol(cm.eqtl.snp.gr$gene), all.genes)
shared.egenes <- intersect(ensembl.to.symbol(shared.eqtl.snp.gr$gene), all.genes)
  

par(mfrow=c(1,2))
boxplot(srna.exp.mat[cm.specific.egenes, "Cardiomyocyte"], 
        srna.exp.mat[shared.egenes, "Cardiomyocyte"], 
        names = c("CM Specific","CM Shared"),
        ylab='CM log2 TP10k',
        main='t-test p: 6e-4')

boxplot(agg.exp.mat[cm.specific.egenes,"agg"], 
        agg.exp.mat[shared.egenes, "agg"],
        names = c("CM Specific","CM Shared"),
        ylab='Bulk log2 TP10k',
        main='t-test p: 0.1')





