source('R/analysis_utils.R')

#hg38
markers <- readRDS('/project2/gca/aselewa/heart_atlas_project/ArchR/ArchR_heart_latest_noAtrium/PeakCalls/DA_MARKERS_FDRP_1_log2FC_1.rds')
cm.markers <- hg38ToHg19(markers$Cardiomyocyte)
seqlevelsStyle(cm.markers) <- "NCBI"

peak.set.hg19 <- hg38ToHg19(peak.set)
seqlevelsStyle(peak.set.hg19) <- "NCBI"

#hg19
finemap.res <- readRDS('GWAS/finemapping/aFib_Finemapped_GeneMapped_ActivePromoter_07222021.gr.rds')
finemap.res <- finemap.res[!duplicated(finemap.res$snp),]
high.pip.snps <- finemap.res[finemap.res$pip>0.5,]
low.pip.snps <- finemap.res[finemap.res$pip<0.01,]


high.pip.snps.gr <- GRanges(seqnames = high.pip.snps$chr, 
                            ranges = IRanges(start = high.pip.snps$pos, end = high.pip.snps$pos),
                            snp = high.pip.snps$snp)

low.pip.snps.gr <- GRanges(seqnames = low.pip.snps$chr, 
                        ranges = IRanges(start = low.pip.snps$pos, end = low.pip.snps$pos),
                        snp = low.pip.snps$snp)

# load all SNPs with MAF > 5%
dbsnp150 <- rtracklayer::import('/project2/xinhe/shared_data/dbsnp150/dbsnp_150_maf05_snpsOnly.vcf.gz') #hg19
dbsnp150.gr <- SummarizedExperiment::rowRanges(dbsnp150)
dbsnp150.gr$SNP_id <- rownames(VariantAnnotation::info(dbsnp150))
seqlevelsStyle(dbsnp150.gr) <- 'UCSC'

#overlap random snps with OCRs
high.pip.snps.gr <- subsetByOverlaps(high.pip.snps.gr, cm.markers)
low.pip.snps.gr <- subsetByOverlaps(low.pip.snps.gr, cm.markers)


nreps <- 15
snp.list <- list()
for(i in 1:nreps){
    snp.list[[i]] <- low.pip.snps.gr[sample(1:length(low.pip.snps.gr), size = length(high.pip.snps.gr), replace = F),]
}

fgt.chip <- readr::read_tsv('ENCODE/FGT_ChIP_lifted_from_mm10.bed', col_names = F)
fgt.chip.gr <- GRanges(seqnames = sub('chr','',fgt.chip$X1), ranges = IRanges(start = fgt.chip$X2, end = fgt.chip$X3), type=fgt.chip$X4)

h3k <- readr::read_tsv('ENCODE/H3k27ac_gwas_hg19/hg19_mapped/H3K27ac_heart_concat.bed', col_names = F)
h3k.gr <- GRanges(seqnames = sub('chr','',h3k$X1), ranges = IRanges(start = h3k$X2, end = h3k$X3))

encode.gr <- list("Fog/Gata4/Tbx5"=fgt.chip.gr, "H3k27ac"=h3k.gr)

high.pip.overlaps <- join_overlap_list(gr.list = encode.gr, X = high.pip.snps.gr)
random.overlaps <- lapply(snp.list, function(x){join_overlap_list(gr.list = encode.gr, X = x)})

high.pip.overlaps.prop <- lapply(high.pip.overlaps, function(x){length(unique(x$snp))/length(high.pip.snps.gr)})

random.overlaps.prop <- unlist(lapply(random.overlaps, function(x){
    sapply(x, function(y){length(unique(y$snp))/length(snp.list[[1]])})
}))

mean.fgt.random.prop <- mean(random.overlaps.prop[names(random.overlaps.prop)=="Fog/Gata4/Tbx5"])
sd.fgt.random.prop <- sd(random.overlaps.prop[names(random.overlaps.prop)=="Fog/Gata4/Tbx5"])/sqrt(length(snp.list[[1]]))

mean.h3.random.prop <- mean(random.overlaps.prop[names(random.overlaps.prop)=="H3k27ac"])
sd.h3.random.prop <- sd(random.overlaps.prop[names(random.overlaps.prop)=="H3k27ac"])/sqrt(length(snp.list[[1]]))

chipseq.df <- data.frame(props = c(high.pip.overlaps.prop$`Fog/Gata4/Tbx5`, high.pip.overlaps.prop$H3k27ac, 
                                   mean.fgt.random.prop, mean.h3.random.prop),
                         type = rep(c("Nkx2-5/Gata4/Tbx5","H3k27ac"), 2),
                         SNPs = rep(c("GWAS SNPs in CM-OCRs (PIP > 0.5)", "GWAS SNPs in CM-OCRs (PIP < 0.01)"), each = 2),
                         sd = c(NA, NA, sd.fgt.random.prop, sd.h3.random.prop))
chipseq.df$SNPs <- factor(chipseq.df$SNPs, levels = c("GWAS SNPs in CM-OCRs (PIP > 0.5)", "GWAS SNPs in CM-OCRs (PIP < 0.01)"))

pdf('manuscript_figures/ChIP_seq_PIP50_overlap.pdf', width=8, height=6)
ggplot(chipseq.df, aes(x=type, y=props, fill=SNPs)) + 
    geom_bar(stat='identity', position='dodge') + 
    ggClean() + 
    ylab('Proportion of SNPs') + 
    xlab('') +
    scale_fill_brewer(palette = 'Set2') + 
    geom_errorbar(aes(ymin=props-(1.96*sd), ymax=props+(1.96*sd)), width=0.1,
                  position=position_dodge(.9)) + coord_cartesian(ylim = c(0, 1))
dev.off()

