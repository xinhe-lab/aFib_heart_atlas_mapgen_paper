library(GenomicRanges)
library(ggplot2)
library(liftOver)
library(tidyverse)
setwd('/project2/gca/aselewa/heart_atlas_project/')
source('R/analysis_utils.R')

satac <- ArchR::loadArchRProject('ArchR/ArchR_heart_latest_noAtrium/')
peak.set <- satac@peakSet

enhancer.coaccess <- readRDS('ArchR/ArchR_heart_latest_noAtrium/CoAccessibility/Coaccess_CMs_enhancer_corr_cut_-1_maxDist_1Mb_hg19.gr.rds') %>% as_tibble() %>% 
    dplyr::select(enhancer_chr, enhancer_start, enhancer_end, correlation, idx, gene_name) 
promoter.coacc <- readRDS('ArchR/ArchR_heart_latest_noAtrium/CoAccessibility/Coaccess_CMs_promoters_corr_cut_-1_maxDist_1Mb_hg19.gr.rds') %>% as_tibble() %>% 
    dplyr::select(promoter_chr, promoter_start, promoter_end, idx)
coacc.df <- enhancer.coaccess %>% inner_join(., promoter.coacc, on = 'idx')

pchic.enhancer.gr <- readRDS('HiC/iPSC_CM_pcHiC_protein_Hg19.gr.rds')
pchic.promoter <- GRanges(seqnames = pchic.enhancer.gr$promoter_chr, ranges = IRanges(start = pchic.enhancer.gr$promoter_start, end = pchic.enhancer.gr$promoter_end),
                           sample = pchic.enhancer.gr$sample)

cor.cutoffs <- list("<0.1"=c(-0.01,0.1), "0.1-0.2"=c(0.1, 0.2), "0.2-0.4"=c(0.2,0.4),"0.4+"=c(0.4,1))
dist.cutoffs <- c(0, seq(100, 500, 100))*1000
res <- list()
overlap.list <- list()

for(i in 1:length(cor.cutoffs)){
    
    curr.cut <- cor.cutoffs[[i]]
    hic.odds <- rep(0, (length(dist.cutoffs)-1))
    overlaps <- rep(0, (length(dist.cutoffs)-1))
    
    coacc.df.high <- coacc.df[(coacc.df$correlation > curr.cut[1]) & (coacc.df$correlation <= curr.cut[2]),]
    dists <- abs(coacc.df.high$enhancer_start-coacc.df.high$promoter_start)
   # rand.dists <- abs(rand.coacc.df$enhancer_start - rand.coacc.df$promoter_start)
    
    for(j in 1:(length(dist.cutoffs)-1)){
        
        coacc.df.high.dist <- coacc.df.high[(dists > dist.cutoffs[j]) & (dists <= dist.cutoffs[j+1]), ]

        enhancer.coacc.gr <- GRanges(seqnames = coacc.df.high.dist$enhancer_chr, 
                                     ranges = IRanges(start = coacc.df.high.dist$enhancer_start, 
                                                      end = coacc.df.high.dist$enhancer_end), 
                                     corr = coacc.df.high.dist$correlation, 
                                     gene = coacc.df.high.dist$gene_name)
        promoter.coacc.gr <- GRanges(seqnames = coacc.df.high.dist$promoter_chr, 
                                     ranges = IRanges(start = coacc.df.high.dist$promoter_start, 
                                                      end = coacc.df.high.dist$promoter_end))
        
        enhancer.hits <- GenomicRanges::findOverlaps(query = enhancer.coacc.gr, subject = pchic.enhancer.gr, maxgap = 1500)
        promoter.hits <- GenomicRanges::findOverlaps(query = promoter.coacc.gr, subject = pchic.promoter, maxgap = 1500)
        
        ehits.str <- paste0(queryHits(enhancer.hits),'-',subjectHits(enhancer.hits))
        phits.str <- paste0(queryHits(promoter.hits),'-',subjectHits(promoter.hits))
        
        ehitsIn <- queryHits(enhancer.hits)
        ehitsIn <- ehitsIn[ehits.str %in% phits.str]
        ehitsIn <- unique(ehitsIn)
        
        propOverlap <- length(ehitsIn)/length(enhancer.coacc.gr)
        overlaps[j] <- propOverlap
        
        hic.odds[j] <- propOverlap
    }
    res[[names(cor.cutoffs)[i]]] <- hic.odds
    overlap.list[[names(cor.cutoffs)[i]]] <- overlaps
}

res.df <- data.frame(overlap = unlist(overlap.list, use.names = F), corr.bin = rep(names(res), lengths(res)), dist.bin = rep(dist.cutoffs[-1], length(cor.cutoffs)))
res.df$corr.bin <- factor(res.df$corr.bin, levels = unique(res.df$corr.bin))

pdf('manuscript_figures/coacc_CMs_overlap_dist_bin_in_iPSC_pcHiC.pdf', width=8, height=5)
ggplot(res.df, aes(x = dist.bin/1000, y = overlap, color=corr.bin)) + 
    geom_point() + 
    geom_line() + 
    ggClean() + 
    xlab('Distance b/w enhancer-promoter') + 
    ylab('Prop. of Coaccessible Peaks in PC HiC') + 
    scale_color_manual(values = c("black","mediumpurple1","royalblue1","skyblue3"))
dev.off()

base.line <- res$`<0`
res.enrich.df <- res.df %>% group_by(corr.bin) %>% mutate(enrich = overlap/base.line) %>% filter(corr.bin != "<0")

pdf('manuscript_figures/coacc_CMs_enrich_dist_bin_in_pcHiC.pdf', width=8, height=5)
ggplot(res.enrich.df, aes(x = dist.bin/1000, y = enrich, color=corr.bin)) + geom_point() + geom_line() + ggClean() + 
    scale_color_brewer(palette = "Blues") + xlab('Distance (kb)') + ylab('Enrichment') +
    geom_hline(yintercept = 1, linetype='dashed', color='grey')
dev.off()

prop.df <- data.frame(overlap = c(propOverlap, randPropOverlap), sample = rep(c(1,2,3), 2), type = c(rep("corr > 0.5", 3), rep("Random Links", 3)))

pdf('manuscript_figures/coacc_CMs_corr10_in_pcHiC.pdf', width=10, height=8)
ggplot(prop.df, aes(x = sample, y = overlap, fill=type)) + geom_bar(stat = 'identity', position = 'dodge') + ylab('Proportion of Co-access Links in PC HiC') +
    scale_fill_brewer(palette = "Set2") + ggClean()
dev.off()
