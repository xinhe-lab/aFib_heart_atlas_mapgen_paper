library(ArchR)
library(clusterProfiler)
library(tidyverse)

require(org.Hs.eg.db)
require(plyranges)
setwd('/project2/gca/aselewa/heart_atlas_project/')
source('R/analysis_utils.R')

# load a bunch of stuff
satac <- suppressMessages(loadArchRProject('ArchR/ArchR_heart_latest_noAtrium/'))

peaks <- getPeakSet(satac)
peaks$peakID <- GRToString(peaks)
peak.markers <- readRDS('ArchR/ArchR_heart_latest_noAtrium/PeakCalls/DA_MARKERS_FDRP_1_log2FC_1.rds')

other.ranges <- readRDS('ArchR/ArchR_heart_latest_noAtrium/CoAccessibility/Coacc_ENHANCERS_AllCellTypes_overlapCutoff50_k200_corr_cut_-1_maxDist_1Mb_hg38.gr.rds')
other.ranges <- other.ranges[other.ranges$correlation > 0.5,]

#annotate all peaks with DA test results
celltype_ideal_order <- c("Cardiomyocyte","Smooth Muscle","Pericyte","Endothelial","Fibroblast","Neuronal", "Lymphoid","Myeloid")
peak.markers <- lapply(celltype_ideal_order, function(x){peak.markers[[x]]})
names(peak.markers) <- celltype_ideal_order
peak.markers.str <- unlist(lapply(peak.markers, function(x){GRToString(x)}), use.names = F)
peaks.celltypes <- data.frame(peakID=peak.markers.str, cellType = factor(rep(names(peak.markers), lengths(peak.markers)), levels = names(peak.markers)))

peak.info.df <- peaks %>% as_tibble() %>% left_join(., peaks.celltypes, on = "peakID")
levels(peak.info.df$cellType) <- c(levels(peak.info.df$cellType), "Shared")
peak.info.df$cellType[is.na(peak.info.df$cellType)] <- "Shared"

# GO
enrich.res <- list()
for(ct in celltype_ideal_order){

  curr <- peak.info.df[peak.info.df$cellType == ct,]
  curr.coacc.genes <- other.ranges$coacc_gene_name[other.ranges$peakID %in% curr$peakID]

  enrich.res[[ct]] <- enrichGO(gene = unique(curr.coacc.genes),
                               OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                               keyType = "SYMBOL",
                               ont = "ALL",
                               qvalueCutoff = 0.05)
}

saveRDS(enrich.res, file = 'misc/enrichGO_Coaccessible_Genes_results.rds')

enrich.res <- readRDS('misc/enrichGO_Coaccessible_Genes_results.rds')
enrich.res <- enrich.res[c(1,4,5,8)]

# Add fold-change results
for (ct in names(enrich.res)){
  enrich_df <- enrich.res[[ct]]@result
  tg_ratio <- enrich_df %>% pull(GeneRatio)
  tg_mat <- do.call(rbind, strsplit(tg_ratio, split = "/"))
  enrich_df$tg_1 <- as.numeric(tg_mat[, 1])
  enrich_df$tg_2 <- as.numeric(tg_mat[, 2])

  bg_ratio <- enrich_df %>% pull(BgRatio)
  bg_mat <- do.call(rbind, strsplit(bg_ratio, split = "/"))
  enrich_df$bg_1 <- as.numeric(bg_mat[, 1])
  enrich_df$bg_2 <- as.numeric(bg_mat[, 2])

  enrich.res[[ct]]@result <- enrich_df %>% mutate(FoldChange = (tg_1/tg_2) / (bg_1/bg_2))
}

# arrange by fold change, pvalue and select relevant columns
top.GO.ids <- list()
for (ct in names(enrich.res)){
  top.GO.ids[[ct]] <- enrich.res[[ct]]@result %>%
    filter(ONTOLOGY == "BP") %>%
    arrange(-FoldChange, pvalue) %>% distinct(geneID, .keep_all = T) %>%
    head(200) %>%
    dplyr::select(ONTOLOGY, ID, Description, FoldChange, pvalue, qvalue, geneID) %>%
    mutate(celltype = ct)
}

for (ct in names(enrich.res)){
  interest.GO.ids <- top.GO.ids[[ct]]$ID
  top.qvalue.mat <- matrix(nrow = length(interest.GO.ids),
                           ncol = length(enrich.res))
  rownames(top.qvalue.mat) <- interest.GO.ids
  colnames(top.qvalue.mat) <- names(enrich.res)
  for (ct2 in names(enrich.res)){
    match.indx <- match(interest.GO.ids, enrich.res[[ct2]]$ID)
    top.qvalue.mat[, ct2] <- enrich.res[[ct2]]$FoldChange[match.indx] #
  }
  top.qvalue.mat[is.na(top.qvalue.mat)] <- 1
  neg.log10.qvalue.df <- reshape2::melt(top.qvalue.mat, value.name = "FC") %>% #
    rename(GO_term = Var1, Cell_type = Var2)
  plot_out <- ggplot(neg.log10.qvalue.df,
                     aes(x = Cell_type, y = GO_term, size = FC, color = Cell_type)) + #
    geom_point() +
    labs(title = paste("Top GO BP terms in", ct)) +
    ggClean()
  print(plot_out)
}

top.GO.df <- do.call(rbind, top.GO.ids) %>% distinct(ID, Description)

# interest.GO.ids <- c("GO:0060047", "GO:0030048", "GO:0007160", "GO:0043542", "GO:0030199", "GO:0010001", "GO:0030098", "GO:0042119") # By pvalue
#interest.GO.ids <- c("GO:0055003", "GO:1905065", "GO:0038166", "GO:0060312", "GO:0030199", "GO:0099560", "GO:0043379", "GO:0002281") # By fold change
#interest.GO.ids <- c("GO:0031033", "GO:0030049","GO:0060837","GO:0060055","GO:1901201","GO:0001941","GO:0002291","GO:0001780")
interest.GO.ids <- c(top.GO.ids$Cardiomyocyte$ID[c(9,11,14,20,21)], top.GO.ids$Endothelial$ID[c(6,13,15,38,82)], top.GO.ids$Fibroblast$ID[c(19,27,29,30,37)],
                     top.GO.ids$Myeloid$ID[c(2,13,19,22,26)])

interest.GO.terms <- top.GO.df$Description[match(interest.GO.ids, top.GO.df$ID)]

top.qvalue.mat <- matrix(nrow = length(interest.GO.ids),
                         ncol = length(enrich.res))
rownames(top.qvalue.mat) <- interest.GO.terms
colnames(top.qvalue.mat) <- names(enrich.res)

top.FC.mat <- top.qvalue.mat

for (ct in names(enrich.res)){
  match.indx <- match(interest.GO.ids, enrich.res[[ct]]$ID)
  top.FC.mat[, ct] <- enrich.res[[ct]]$FoldChange[match.indx] #
  top.qvalue.mat[, ct] <- -log10(enrich.res[[ct]]$qvalue[match.indx])
}

top.qvalue.mat[is.na(top.qvalue.mat)] <- 0
neg.log10.qvalue.df <- reshape2::melt(top.qvalue.mat, value.name = "mlogqval") %>% #
  rename(GO_term = Var1, Cell_type = Var2)

top.FC.mat[is.na(top.FC.mat)] <- 0.1
qval.foldchange.df <- reshape2::melt(top.FC.mat, value.name = "FC") %>% #
  rename(GO_term = Var1, Cell_type = Var2) %>% left_join(., neg.log10.qvalue.df) %>%
  mutate(GO_term = factor(GO_term, levels = rev(unique(interest.GO.terms))))


pdf('manuscript_figures/figure3/Fig3_EnrichGO_Coaccess.pdf', width=12, height=9)

ggplot(qval.foldchange.df) +
  geom_point(aes(x = Cell_type, y = GO_term, size = FC, color = mlogqval)) +
  scale_size_continuous(range = c(1, 10)) +
  ggClean(rotate_axis = T) +
  scale_color_gradientn(colors = c('navy',"yellow",'red')) +
  geom_hline(yintercept = c(5.5,10.5,15.5), linetype='dashed')

dev.off()



# enrichGO analysis of high pip genes

gene.map <- read_csv('GWAS/aFib_Finemapped_GenePIP_0.1_ActiveProm_07222021.csv')
high.pip.genes <- unique(gene.map$gene_name[gene.map$gene_pip > 0.5])

gene.map.enrich <- enrichGO(gene = high.pip.genes,
                            OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                            keyType = "SYMBOL",
                            ont = "ALL",
                            qvalueCutoff = 0.05)

enrich_df <- gene.map.enrich@result
tg_ratio <- enrich_df %>% pull(GeneRatio)
tg_mat <- do.call(rbind, strsplit(tg_ratio, split = "/"))
enrich_df$tg_1 <- as.numeric(tg_mat[, 1])
enrich_df$tg_2 <- as.numeric(tg_mat[, 2])

bg_ratio <- enrich_df %>% pull(BgRatio)
bg_mat <- do.call(rbind, strsplit(bg_ratio, split = "/"))
enrich_df$bg_1 <- as.numeric(bg_mat[, 1])
enrich_df$bg_2 <- as.numeric(bg_mat[, 2])

gene.map.enrich@result <- enrich_df %>% mutate(FoldChange = (tg_1/tg_2) / (bg_1/bg_2))

# saveRDS(gene.map.enrich, file = 'misc/gene_mapping_enrich_GO.rds')

go.res <- gene.map.enrich@result %>% arrange(-FoldChange) %>% distinct(geneID, .keep_all = T)
top.bp.go <- go.res %>% filter(ONTOLOGY == "BP") %>% slice(1:10)
top.mf.go <- go.res %>% filter(ONTOLOGY == "MF") %>% slice(1:10)

top.bp.go$Description <- factor(top.bp.go$Description, levels = rev(top.bp.go$Description))
top.mf.go$Description <- factor(top.mf.go$Description, levels = rev(top.mf.go$Description))

# pdf('manuscript_figures/GWAS_bp_genes_GOenrich.pdf',width=12,height=8)

ggplot(top.bp.go, aes(x = FoldChange, y = Description)) +
  geom_bar(stat='identity', fill="grey") +
  scale_size_continuous(range = c(1, 10)) +
  ggClean(rotate_axis = T) +
  xlab('Enrichment')

# dev.off()

# pdf('manuscript_figures/GWAS_mf_genes_GOenrich.pdf',width=16,height=8)

ggplot(top.mf.go, aes(x = FoldChange, y = Description)) +
  geom_bar(stat='identity', fill = "slategray1") +
  scale_size_continuous(range = c(1, 10)) +
  ggClean(rotate_axis = T) +
  xlab('Enrichment')

# dev.off()

# go.res %>% arrange(-FoldChange) %>% write_csv('High_PIP_Genes_GO_Results.csv')

