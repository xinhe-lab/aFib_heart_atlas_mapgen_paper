library(gkmSVM)
library(tidyverse)
library(liftOver)
library(BSgenome.Hsapiens.UCSC.hg19.masked)
setwd('/project2/gca/aselewa/heart_atlas_project/')

markers <- readRDS('ArchR/ArchR_heart_latest_noAtrium/PeakCalls/DA_MARKERS_FDRP_1_log2FC_1.rds')

path <- system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch <- import.chain(path)
markers.hg19 <- lapply(markers, function(x){unlist(liftOver(x, ch))})
for(i in 1:length(markers.hg19)){
    markers.hg19[[i]] <- markers.hg19[[i]][width(markers.hg19[[i]]) >= 100,]
    seqlevelsStyle(markers.hg19[[i]]) <- "UCSC"
}
system('mkdir -p GWAS/bed_annotations_hg19')
lapply(names(markers.hg19), function(x){rtracklayer::export(markers.hg19[[x]], 
                                                            format = 'bed', 
                                                            con = paste0('GWAS/bed_annotations_hg19/',x,'_narrowPeaks_UCSC.bed'))})

system('mkdir -p GWAS/deltaSVM')

annots <- list.files('GWAS/bed_annotations_hg19', 'UCSC.bed', full.names = T)
for(b in annots){
    fname <- sub(' ', '', sub('_.*', '', basename(b)))
    genNullSeqs(inputBedFN = b, 
                outputPosFastaFN = paste0('GWAS/deltaSVM/',fname,'_posSet.fa'), 
                outputNegFastaFN = paste0('GWAS/deltaSVM/',fname,'_negSet.fa'))
}

system('/project2/gca/software/miniconda3/envs/ldsc/bin/python /project2/gca/software/lsgkm/scripts/nrkmers.py 11 GWAS/deltaSVM/n11mers.fa')

system('/project2/gca/software/lsgkm/src/gkmtrain GWAS/deltaSVM/Cardiomyocyte_posSet.fa GWAS/deltaSVM/Cardiomyocyte_negSet.fa GWAS/deltaSVM/CM -T 16 -m 100000')

system('/project2/gca/software/lsgkm/src/gkmpredict GWAS/deltaSVM/n11mers.fa GWAS/deltaSVM/CM.model.txt GWAS/deltaSVM/CM_n11mers_weights.txt -T 4')

# Prepare reference and alt allele fasta seqs
finemap.res <- readRDS('GWAS/finemapping/aFib_Finemapped.tble.rds')
high.pip.snps <- finemap.res[finemap.res$susie_pip > 0.2, ]

alt.fa <- rep("", nrow(high.pip.snps))
ref.fa <- rep("", nrow(high.pip.snps))
offset <- 9
for(i in 1:nrow(high.pip.snps)){
    chrom <- paste0("chr",high.pip.snps$chr[i])
    pos <- high.pip.snps$pos[i]
    ref <- high.pip.snps$a0[i]
    alt <- high.pip.snps$a1[i]
    
    var_id <- paste0(">",chrom,'_',pos,'_',ref,'_',alt)
    left_flank <- getSeq(BSgenome.Hsapiens.UCSC.hg19.masked, chrom, pos-offset, pos-1)
    right_flank <- getSeq(BSgenome.Hsapiens.UCSC.hg19.masked, chrom, pos+1, pos+offset)
    
    ref.fa[i] <- paste0(var_id,'\n',left_flank,ref,right_flank)
    alt.fa[i] <- paste0(var_id,'\n',left_flank,alt,right_flank)
}

writeLines(ref.fa, con = 'GWAS/deltaSVM/ref.fa')
writeLines(alt.fa, con = 'GWAS/deltaSVM/alt.fa')

# Run model to predict effects
system('perl /project2/gca/software/deltasvm_script/deltasvm.pl GWAS/deltaSVM/ref.fa GWAS/deltaSVM/alt.fa GWAS/deltaSVM/CM_n11mers_weights.txt GWAS/deltaSVM/PIP20_CM_ATAC_pred.txt')

snp.delta <- read_tsv('GWAS/deltaSVM/PIP20_CM_ATAC_pred.txt', col_names = F)

# Run model on permuted weights to get z-scores
system('mkdir -p GWAS/deltaSVM/permutations')

kmer.weights <- read_tsv('GWAS/deltaSVM/CM_n11mers_weights.txt', col_names = F)
N_permutes <- 100
for(i in 1:N_permutes){
    transform(kmer.weights, X2 = sample(kmer.weights$X2)) %>% write_tsv(paste0('GWAS/deltaSVM/permutations/CM_n11mers_weights_permut_',i,'.txt'), col_names = F)
    system(paste0('perl /project2/gca/software/deltasvm_script/deltasvm.pl GWAS/deltaSVM/ref.fa GWAS/deltaSVM/alt.fa GWAS/deltaSVM/permutations/CM_n11mers_weights_permut_',i,'.txt GWAS/deltaSVM/permutations/PIP20_CM_ATAC_pred_permute_',i,'.txt'))
}

# analyze the deltaSVM z-scores
svm.scores <- list()
files <- list.files('GWAS/deltaSVM/permutations/', 'PIP20_', full.names = T)
for(f in 1:length(files)){
    svm.scores[[f]] <- read_tsv(files[f], col_names = F)$X2
}
svm.null.mat <- Reduce(cbind, svm.scores)
snp.delta$zscore <- ((snp.delta$X2)-rowMeans(svm.null.mat))/matrixStats::rowSds(svm.null.mat)
snp.delta$pvals <- p.adjust(2*(1-pnorm(abs(snp.delta$zscore))), method = "BH")

snp.annots <- finemappeR::annotator(gwas = high.pip.snps, annotations = list.files('GWAS/bed_annotations_hg19/', 'Cardiomyocyte_narrowPeaks.bed', full.names = T))
snp.annots$deltaSVM_zscore <- snp.delta$zscore
snp.annots$padjust <- snp.delta$pvals

# run deltaSVM on some TFs
system('/project2/gca/software/lsgkm/src/gkmpredict GWAS/deltaSVM/n11mers.fa /project2/gca/software/deltaSVM/gkmsvm_models/ESRRG_eDBD_4_KS_TACCTT40NCGA.subsample_20000.10mer.model.txt GWAS/deltaSVM/CM_n11mers_ESRRG_weights.txt -T 4')
system('perl /project2/gca/software/deltasvm_script/deltasvm.pl GWAS/deltaSVM/ref.fa GWAS/deltaSVM/alt.fa GWAS/deltaSVM/CM_n11mers_ESRRG_weights.txt GWAS/deltaSVM/PIP20_CM_ESRRG_pred.txt')

esrrg.delta <- read_tsv('GWAS/deltaSVM/PIP20_CM_ESRRG_pred.txt', col_names = F)

plot(esrrg.delta$X2, snp.delta$X2)



