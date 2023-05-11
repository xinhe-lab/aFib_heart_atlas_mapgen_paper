source('R/analysis_utils.R')
library(motifbreakR)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)

finemap.res <- readRDS('GWAS/finemapping/aFib_Finemapped_GeneMapped.tble.rds')
finemap.res <- finemap.res[!duplicated(finemap.res$snp),]
high.pip.snps <- finemap.res[finemap.res$pip>0.5,]
high.pip.snps <- high.pip.snps[high.pip.snps$annots!="Other" & high.pip.snps$annots!="Non-CM ATAC",]$snp

markers <- readRDS('/project2/gca/aselewa/heart_atlas_project/ArchR/ArchR_heart_latest_noAtrium/PeakCalls/DA_MARKERS_FDRP_1_log2FC_1.rds')

# load all SNPs with MAF > 5%
dbsnp150 <- rtracklayer::import('/project2/xinhe/shared_data/dbsnp150/dbsnp_150_maf05_snpsOnly.vcf.gz') #hg19
dbsnp150.gr <- SummarizedExperiment::rowRanges(dbsnp150)
dbsnp150.gr$SNP_id <- rownames(VariantAnnotation::info(dbsnp150))
seqlevelsStyle(dbsnp150.gr) <- 'UCSC'

#overlap random snps with OCRs
dbsnp150.in.OCR <- subsetByOverlaps(dbsnp150.gr, markers$Cardiomyocyte)

nreps <- 15
snp.list <- list()
for(i in 1:nreps){
    snp.list[[i]] <- sample(dbsnp150.in.OCR$SNP_id, size = length(high.pip.snps), replace = F)
}

# start motifbreakR
cardiac.tfs <- c("TBX5","GATA4","PITX2","ETV1", "NKX2-5","PRRX1","MEF2A","MEF2C","TEAD1","NR2F2","HAND1")

res <- query(MotifDb, andStrings = c("hsapiens"), orStrings = cardiac.tfs)
res <- res[!duplicated(mcols(res)$geneSymbol)]
res <- res[mcols(res)$geneSymbol %in% cardiac.tfs,]

hand2.pwm <- read.table('misc/HAND2_motif.txt', sep = '')
colnames(hand2.pwm) <- seq(1:10)
rownames(hand2.pwm) <- c("A","C","G","T")
hand2.pwm <- t(t(hand2.pwm)/colSums(hand2.pwm))
res$HAND2 <- hand2.pwm
mcols(res)["HAND2",] <- "HAND2"

rand.motif.breaks.res <- list()
for(i in 1:nreps){
    variants <- snps.from.rsid(rsid = snp.list[[i]], 
                               dbSNP = SNPlocs.Hsapiens.dbSNP142.GRCh37::SNPlocs.Hsapiens.dbSNP142.GRCh37,
                               search.genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
    variants <- variants[!duplicated(variants$SNP_id),]
    print(length(variants))
    rand.motif.breaks.res[[i]] <- motifbreakR(snpList = variants, 
                                show.neutral = TRUE,
                                filterp = TRUE,
                                pwmList = res, 
                                threshold = 1e-3,
                                method = "ic",
                                bkg = c(A=0.25, C=0.25, G=0.25, T=0.25), 
                                BPPARAM = BiocParallel::MulticoreParam(workers = 16))    
}

saveRDS(rand.motif.breaks.res, 'misc/Linsin_tfs_motifbreakR_randomSNPs_in_OCRs_nrep15_PIP50.gr.rds')

# run high pip case
variants <- snps.from.rsid(rsid = high.pip.snps,
                           dbSNP = SNPlocs.Hsapiens.dbSNP142.GRCh37::SNPlocs.Hsapiens.dbSNP142.GRCh37,
                           search.genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
drops <- c("rs3731746:T","rs60632610:A","rs60632610:G")
variants <- variants[!(names(variants) %in% drops),]
motif.breaks.high.pip <- motifbreakR(snpList = variants, 
                            show.neutral = TRUE,
                            filterp = TRUE,
                            pwmList = res, 
                            threshold = 1e-3,
                            method = "ic",
                            bkg = c(A=0.25, C=0.25, G=0.25, T=0.25), 
                            BPPARAM = BiocParallel::MulticoreParam(workers = 16))

saveRDS(motif.breaks.high.pip, 'misc/Linsin_tfs_motifbreakR_PIP50.gr.rds')

# number of snps in a motif

rand.snp.props <- lapply(rand.motif.breaks.res, function(x){length(unique(x$SNP_id))/length(high.pip.snps)})
z <- unlist(rand.snp.props)
z.mean <- mean(z)
z.sd <- sd(z)/nreps

nsnp.df <- data.frame(prop.snp = c(length(unique(motif.breaks.high.pip$SNP_id))/length(snp.list[[1]]),
                                   z.mean),
                      SNPs = c("GWAS SNPs \n(PIP>0.5)","Random \nCM-OCR SNPs"),
                      sd = c(NA, z.sd))

pdf('SNP_props_TF_motifbreakR_randomOCR_PIP50.pdf', width = 8, height=6)
ggplot(nsnp.df, aes(x=SNPs, y=prop.snp, fill=SNPs)) +  geom_bar(stat='identity') + ggClean() + LegendOff() + scale_fill_brewer(palette = 'Set2') +
    geom_errorbar(aes(ymin=prop.snp-(1.96*sd), ymax=prop.snp+(1.96*sd)), width=.1,
                  position=position_dodge(.9)) +
    ylab('Prop SNPs in Motif')
dev.off()

# bin and plot allele diff score
bks <- c(0.001, 1.5, 2) # bins based on quantiles of allele allele diff scores. 
labs <- c("weak","strong")
bin.count.rand <- lapply(rand.motif.breaks.res, function(x){table(cut(abs(x$alleleDiff), breaks = bks, labels = labs))})
z <- unlist(bin.count.rand)
weak.mean <- mean(z[names(z)=="weak"])
weak.std <- sd(z[names(z)=="weak"])/sqrt(nreps)

strong.mean <-  mean(z[names(z)=="strong"])
strong.std <-  sd(z[names(z)=="strong"])/sqrt(nreps)

bin.count.high.pip <- table(cut(abs(motif.breaks.high.pip$alleleDiff), breaks = bks, labels = labs))

motifbreak.compare <- data.frame(counts=c(bin.count.high.pip, c(weak.mean, strong.mean)),
                                 sd=c(NA, NA, weak.std, strong.std),
                                 SNPs = rep(c("GWAS SNPs \n(PIP>0.5)","Random \nCM-OCR SNPs"), each = 2),
                                 bin = factor(rep(labs,2), levels = labs))

pdf('TF_motifbreakR_randomOCR_PIP50.pdf', width = 8, height=8)
ggplot(motifbreak.compare, aes(x=bin, y=counts, fill=SNPs)) + geom_bar(stat='identity', position = 'dodge') + ggClean() + ylab('Number of SNP-Motif Pairs') +
    xlab('Binned Motif Disrupt Score') +
    geom_errorbar(aes(ymin=counts-(1.96*sd), ymax=counts+(1.96*sd)), width=.1,
                  position=position_dodge(.9)) + 
    scale_fill_brewer(palette = 'Set2')
dev.off()

