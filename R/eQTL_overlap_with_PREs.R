library(ArchR)
library(tidyverse)
setwd('/project2/gca/aselewa/heart_atlas_project/')
source('R/analysis_utils.R')
palette <- readRDS('notebooks/palette.rds')

satac <- suppressMessages(loadArchRProject('ArchR/ArchR_heart_latest_noAtrium/'))
peak.markers <- readRDS('ArchR/ArchR_heart_latest_noAtrium/PeakCalls/DA_MARKERS_FDRP_1_log2FC_1.rds')
celltype_ideal_order <- c("Cardiomyocyte","Smooth Muscle","Pericyte","Endothelial","Fibroblast","Neuronal", "Lymphoid","Myeloid")
peak.set <- satac@peakSet
peak.set$peakID <- GRToString(peak.set)

eqtl.snp.gr <- readRDS('eQTL_enrich/gtex_finemapping/HighPIP_eQTLs_w_LFSR.gr.rds')

rand.snps <- readLines('matched_SNPs/eQTL_top_pip_hg19_5batches/snpsnap_match_hg38.txt') # load random SNPs
rand.gr <- StringToGR(rand.snps)
rand.gr <- rand.gr[seqnames(rand.gr) %in% paste0("chr",1:22),]
seqlevels(rand.gr) <-  paste0("chr",1:22)
rand.gr$snp <- paste0(seqnames(rand.gr),'_',start(rand.gr))

peak.markers.union <- split(peak.set, peak.set$peakType)
eqtl.overlap.peaks <- lapply(peak.markers.union, function(x){plyranges::join_overlap_inner(eqtl.snp.gr, x)})
eqtl.overlap.peaks.prop <- sapply(eqtl.overlap.peaks, function(x){length(x)/length(eqtl.snp.gr)})

random.overlap.peaks <- lapply(peak.markers.union, function(x){plyranges::join_overlap_inner(rand.gr, x)})
random.overlap.peaks.prop <- sapply(random.overlap.peaks, function(x){length(x)/length(rand.gr)})

overlap.peaks.df <- data.frame(overlap.peaks.prop = c(eqtl.overlap.peaks.prop, random.overlap.peaks.prop))
overlap.peaks.df$celltypes <- factor(rep(names(peak.markers.union),2))
#, levels = c("union",celltype_ideal_order))
overlap.peaks.df$snpType <- rep(c("eQTLs","Random"), each = length(peak.markers.union))

pdf('manuscript_figures/figure4/Fig4B_union_peaktype_overlap.pdf', width = 6, height=8)
ggplot(overlap.peaks.df, aes(x = celltypes, y = overlap.peaks.prop, fill = snpType)) + 
    geom_bar(stat='identity', position='dodge') + ggClean(rotate_axis = T) +
    scale_fill_manual(values=c("black","grey")) + 
    xlab('Cell Types') + ylab('% eQTLs')
dev.off()

peak.markers.union <- peak.markers
eqtl.overlap.peaks <- lapply(peak.markers.union, function(x){plyranges::join_overlap_inner(eqtl.snp.gr, x)})
eqtl.overlap.peaks.prop <- sapply(eqtl.overlap.peaks, function(x){length(x)/length(eqtl.snp.gr)})

random.overlap.peaks <- lapply(peak.markers.union, function(x){plyranges::join_overlap_inner(rand.gr, x)})
random.overlap.peaks.prop <- sapply(random.overlap.peaks, function(x){length(x)/length(rand.gr)})

overlap.peaks.df <- data.frame(overlap.peaks.prop = c(eqtl.overlap.peaks.prop, random.overlap.peaks.prop))
overlap.peaks.df$celltypes <- factor(rep(names(peak.markers.union),2), levels = c("union",celltype_ideal_order))
overlap.peaks.df$snpType <- rep(c("eQTLs","Random"), each = length(peak.markers.union))

pdf('manuscript_figures/figure4/Fig4B_celltype_overlap.pdf', width = 10, height=8)

ggplot(overlap.peaks.df, aes(x = celltypes, y = overlap.peaks.prop, fill = snpType)) + 
    geom_bar(stat='identity', position='dodge') + ggClean(rotate_axis = T) +
    scale_fill_manual(values=c("black","grey")) + 
    xlab('Cell Types') + ylab('% eQTLs')

dev.off()