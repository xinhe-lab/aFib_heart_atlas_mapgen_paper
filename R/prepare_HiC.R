library(tidyverse)

gene.annots <- readRDS('genomic_annotations/hg38_gtf_genomic_annots.gr.rds')

# CM HiC data
cm.hic <- read_tsv('HiC/pcHiC_CMs_consensus_set.txt', col_names = F) #hg19

bait.gene.names <- sub(pattern = "\\*.*", replacement = "" , cm.hic$X8)
cm.hic.enhancer.gr <- GRanges(seqnames = cm.hic$X4, IRanges(start = cm.hic$X5, end = cm.hic$X6)) 
cm.hic.enhancer.gr$promoter_chr <- cm.hic$X1
cm.hic.enhancer.gr$promoter_start <- cm.hic$X2
cm.hic.enhancer.gr$promoter_end <- cm.hic$X3
cm.hic.enhancer.gr$gene_name <- bait.gene.names
cm.hic.enhancer.gr$score <- cm.hic$X7

cm.hic.enhancer.protein.gr <- cm.hic.enhancer.gr[cm.hic.enhancer.gr$gene_name %in% gene.annots$genes$gene_name,] # restrict to protein coding genes

saveRDS(cm.hic.enhancer.protein.gr, 'HiC/iPSC_CM_pcHiC_protein_Hg19.gr.rds')
