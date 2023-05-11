library(tidyverse)
library(ArchR)
setwd('/project2/gca/aselewa/heart_atlas_project/')
#fine-mapping results
final.mat <- readRDS('GWAS/finemapping/aFib_Finemapped_GeneMapped_06262021.tble.rds')
high.conf.snp.df <- final.mat %>% 
    dplyr::filter(pip > 0.5) %>% 
    group_by(snp) %>% 
    arrange(-gene_pip) %>% 
    dplyr::slice(1) 

peak.set <- satac@peakSet
peak.set.hg19 <- hg38ToHg19(peak.set)

# random snp source
dbsnp150 <- rtracklayer::import('/project2/xinhe/shared_data/dbsnp150/dbsnp_150_maf05_snpsOnly.vcf.gz') #hg19
dbsnp150.gr <- SummarizedExperiment::rowRanges(dbsnp150)
dbsnp150.gr$SNP_id <- rownames(VariantAnnotation::info(dbsnp150))
seqlevelsStyle(dbsnp150.gr) <- 'UCSC'

#overlap random snps with OCRs
dbsnp150.in.OCR <- subsetByOverlaps(dbsnp150.gr, peak.set.hg19)
seqlevelsStyle(dbsnp150.in.OCR) <- 'UCSC'

nreps <- 50
snp.list <- list()
for(i in 1:nreps){
    snp.list[[i]] <- dbsnp150.in.OCR[sample(1:length(dbsnp150.in.OCR), size = nrow(high.conf.snp.df), replace = F),]
}

markers <- readRDS('ArchR/ArchR_heart_latest_noAtrium/PeakCalls/DA_MARKERS_FDRP_1_log2FC_1.rds')

CM.calls <- hg38ToHg19(readRDS('ArchR/ArchR_heart_latest_noAtrium/PeakCalls/Cardiomyocyte-reproduciblePeaks.gr.rds'))
endo.calls <- hg38ToHg19(readRDS('ArchR/ArchR_heart_latest_noAtrium/PeakCalls/Endothelial-reproduciblePeaks.gr.rds'))
fibro.calls <- hg38ToHg19(readRDS('ArchR/ArchR_heart_latest_noAtrium/PeakCalls/Fibroblast-reproduciblePeaks.gr.rds'))

CMSpecific <- rtracklayer::import('GWAS/annotations_for_finemapping_hg19/CM_specific_peaks_hg19.bed')
seqlevelsStyle(CMSpecific) <- "UCSC"
CMShared <- rtracklayer::import('GWAS/annotations_for_finemapping_hg19/CM_shared_peaks_hg19.bed')
seqlevelsStyle(CMShared) <- "UCSC"
NonCM <- rtracklayer::import('GWAS/annotations_for_finemapping_hg19/non_CM_peaks_hg19.bed')
seqlevelsStyle(NonCM) <- "UCSC"

annots.list <- list(CMSpecific=CMSpecific, CMShared=CMShared, NonCM=NonCM)

genomic.annots <- readRDS('genomic_annotations/hg19_gtf_genomic_annots.gr.rds')
annots.list$Exons <- genomic.annots$exons
annots.list$Promoters <- genomic.annots$promoters
annots.list$UTR <- genomic.annots$UTRs


snp.gr <- GRanges(seqnames = paste0("chr",high.conf.snp.df$chr),
                  ranges = IRanges(start = high.conf.snp.df$pos, end = high.conf.snp.df$pos), 
                  snp = high.conf.snp.df$snp)

snps.overlap.annots <- join_overlap_list(gr.list = annots.list, X = snp.gr)
high.pip.overlaps.prop <- unlist(lapply(snps.overlap.annots, function(x){length(unique(x$snp))/length(snp.gr)}))

random.overlaps <- lapply(snp.list, function(x){join_overlap_list(gr.list = annots.list, X = x)})
random.overlaps.prop <- unlist(lapply(random.overlaps, function(x){
    sapply(x, function(y){length(unique(y$SNP_id))/length(snp.list[[1]])})
}))

mean.vec <- rep(0, length(annots.list))
sd.vec <- c(0, length(annots.list))
annots.names <- names(annots.list)
for(i in 1:length(annots.list)){
    mean.vec[i] <- mean(random.overlaps.prop[names(random.overlaps.prop)==annots.names[i]])
    sd.vec[i] <- sd(random.overlaps.prop[names(random.overlaps.prop)==annots.names[i]])/sqrt(length(snp.list[[1]]))
}

snp.loc.enrich <- data.frame(props = c(high.pip.overlaps.prop, mean.vec),
                         type = rep(annots.names, 2),
                         SNPs = rep(c("GWAS SNPs (PIP > 0.5)", "Random OCR SNPs"), each = length(annots.names)),
                         sd = c(rep(NA, length(annots.names)), sd.vec))
snp.loc.enrich$type <- factor(snp.loc.enrich$type, levels = c("CMSpecific","CMShared","NonCM","Exons","Promoters","UTR"))

ggplot(snp.loc.enrich, aes(x=type, y=props, fill=SNPs)) + geom_bar(stat='identity', position='dodge') + ggClean() + ylab('Proportion of SNPs') + xlab('') +
    scale_fill_brewer(palette = 'Set2') + 
    geom_errorbar(aes(ymin=props-(1.96*sd), ymax=props+(1.96*sd)), width=0.1,
                  position=position_dodge(.9))
