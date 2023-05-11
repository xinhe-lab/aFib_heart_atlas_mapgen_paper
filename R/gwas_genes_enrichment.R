library(tidyverse)

setwd('/project2/gca/aselewa/heart_atlas_project/')
source('R/analysis_utils.R')

srna <- readRDS('seurat/Heart_RNA_Processed_Combined_NoAtrium.rds')
srna.exp <- Seurat::AverageExpression(srna)
avg.exp.mat <- log2(srna.exp$RNA + 1)

gene.annots <- readRDS('genomic_annotations/hg19_gtf_genomic_annots.gr.rds')
genes <- gene.annots$genes

diff.exp.res <- readRDS('seurat/diff_expr_markers.df.rds')
finemap.res <- readRDS('GWAS/aFib_activepromoter.OCRs_enhancerloop.ABC.pcHiC.nearby20kb_dist.50000_genemapping_FinalMat_08052021.rds')

high.pip.snps <- finemap.res %>% distinct(snp, .keep_all=TRUE) %>% group_by(locus) %>% arrange(-pip) %>% dplyr::slice(1)
high.pip.snp.gr <- GRanges(seqnames = paste0("chr",high.pip.snps$chr), ranges = IRanges(start = high.pip.snps$pos, end = high.pip.snps$pos))
gene.view <- read_csv('GWAS/aFib_Finemapped_GeneView_GenePIP_0.1_08062021.csv')
top.pip.genes <- gene.view$gene_name[gene.view$gene_pip > 0.8]

# enrichment of GWAS genes among diff exp genes
hits <- findOverlaps(high.pip.snp.gr, genes, maxgap = 250000)
control.genes <- unique(genes$gene_name[subjectHits(hits)])

cell.types <- unique(diff.exp.res$cluster)
top.pip.genes.overlap <- sapply(cell.types, function(x){
  length(intersect(diff.exp.res$gene[diff.exp.res$cluster==x], top.pip.genes))/length(top.pip.genes)})

control.genes.overlap <- sapply(cell.types, function(x){
  length(intersect(diff.exp.res$gene[diff.exp.res$cluster==x], control.genes))/length(control.genes)})

compare.overlap.df <- data.frame(prop = c(top.pip.genes.overlap, control.genes.overlap), 
                                 celltype = as.character(cell.types),
                                 gene_type = rep(c("PIP>0.5","Control"),each=length(cell.types)))

compare.overlap.df <- compare.overlap.df[compare.overlap.df$celltype %in% c("Cardiomyocyte","Endothelial","Fibroblast"),]

pdf('manuscript_figures/GWAS_Genes_Control_PIP80.pdf',width=8, height=6)
ggplot(compare.overlap.df, aes(x=celltype, y=prop, fill=gene_type)) + geom_bar(stat='identity', position='dodge', width = 0.7) + ggClean(rotate_axis = T) +
  xlab('') + ylab('Prop. Genes Differentially Expressed') + scale_fill_brewer(palette = "Paired")
dev.off()

# control vs. high pip TP10k 
avg.exp.mat.top.pip <- avg.exp.mat[intersect(rownames(avg.exp.mat), top.pip.genes),]
avg.exp.mat.top.pip <- pivot_longer(data = avg.exp.mat.top.pip, cols = everything())
avg.exp.mat.top.pip$gene_type <- "PIP>0.8"

avg.exp.mat.control <- avg.exp.mat[intersect(rownames(avg.exp.mat), control.genes),]
avg.exp.mat.control <- pivot_longer(data = avg.exp.mat.control, cols = everything())
avg.exp.mat.control$gene_type <- "Control"

exp.compare.df <- bind_rows(avg.exp.mat.top.pip,avg.exp.mat.control)
exp.compare.df$value[exp.compare.df$value > 3] <- 3

exp.compare.df <- exp.compare.df[exp.compare.df$name == "Cardiomyocyte",]

pdf('manuscript_figures/GWAS_CM_genes_vs_control_tpm_PIP80.pdf', width=4, height=5)
ggplot(exp.compare.df, aes(x=gene_type, y=value, fill=gene_type)) + geom_boxplot(outlier.size = 0.5) + ggClean(rotate_axis = T) + 
  xlab('') + ylab('log2 TP10k') + LegendOff() + scale_fill_brewer(palette = "Paired")
dev.off()


# enrichment in OMIM genes
all.genes.pips <- finemap.res %>% group_by(gene_name) %>% slice(1)

bks <- c(0,0.1,0.5,0.8,10)
labs <- c("<10%","11%-50%","50%-80%",">80%")
bin.count <- cut(all.genes.pips$gene_pip, breaks = bks, labels = labs)
all.genes.pips$gene_pip_bin <- bin.count

omim.genes <- read_tsv('GWAS/OMIM_GENES_cardiac_arryhthmia_OR_conduction_defect.txt', col_names = F)
colnames(omim.genes) <- c("code","description","gene")

omim.overlap <- all.genes.pips %>% 
  group_by(gene_pip_bin) %>% 
  summarise(prop_omim = length(intersect(gene_name, omim.genes$gene))/dplyr::n())

pdf('manuscript_figures/OMIM_genes_proportion_PIP80.pdf', width=5, height=5)
ggplot(omim.overlap, aes(x=gene_pip_bin, y=prop_omim*100)) + 
  geom_bar(stat='identity', fill='skyblue') + 
  ggClean() + 
  ylab('% OMIM') +
  xlab('Gene PIP Bin')
dev.off()

