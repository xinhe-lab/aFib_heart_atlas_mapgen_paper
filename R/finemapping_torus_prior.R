
## finemapping AF loci with TORUS prior using new ATAC peaks
library(mapgen)
library(dplyr)
suppressMessages(library(GenomicRanges))
suppressMessages(library(data.table))
suppressMessages(library(liftOver))
library(susieR)

setwd('/project2/xinhe/shared_data/aFib_gene_finemapping/pipeline/')

gwas.sumstats <- readRDS(paste0('GWAS/ebi-a-GCST006414_aFib.df.rds'))
head(gwas.sumstats)

# Prepare ATAC peak annotation files -------------------------------------------------------------------------------------------

# load peaks (in hg38)
Categorized_OCRs <- readRDS("/project2/gca/aselewa/heart_atlas_project/ArchR/ArchR_heart_latest_noAtrium/PeakCalls/Categorized_OCRs_all_peaks.gr.rds")
CM.pseudobulk.peaks <- readRDS("/project2/gca/aselewa/heart_atlas_project/ArchR/ArchR_heart_latest_noAtrium/PeakCalls/Cardiomyocyte-reproduciblePeaks.gr.rds")

# CM-specific: significantly differentially upregulated peaks.
CM.specific.peaks <- Categorized_OCRs[Categorized_OCRs$category == "Cardiomyocyte"]
cat(length(CM.specific.peaks), "CM specific peaks.\n")

# CM Non-DA peaks: overlap non-DA peaks with pseudo-bulk peaks
nonDA.peaks <- Categorized_OCRs[Categorized_OCRs$category == "Non-DA"]
CM.nonDA.peaks <- subsetByOverlaps(nonDA.peaks, CM.pseudobulk.peaks)
cat(length(CM.nonDA.peaks), "CM non-DA peaks.\n")

# CM-shared:
shared.peaks <- Categorized_OCRs[Categorized_OCRs$category %in% c("2-3", "4+")]
CM.shared.peaks <- subsetByOverlaps(shared.peaks, CM.pseudobulk.peaks)
cat(length(CM.shared.peaks), "CM shared peaks.\n")
table(CM.shared.peaks$category)

# Non-CM: Anything left in the union set that doesn’t overlap with the above.
CM.idx <- which(countOverlaps(Categorized_OCRs, c(CM.specific.peaks, CM.nonDA.peaks, CM.shared.peaks)) > 0)
other.peaks <- unique(Categorized_OCRs[-CM.idx])
cat(length(other.peaks), "other peaks.\n")

# setequal(c(CM.specific.peaks, CM.nonDA.peaks, CM.shared.peaks, other.peaks), Categorized_OCRs)
# length(unique(c(CM.specific.peaks, CM.nonDA.peaks, CM.shared.peaks, other.peaks)))

peak.sets <- list(`CM_specific_peaks` = CM.specific.peaks,
                  `CM_nonDA_peaks` = CM.nonDA.peaks,
                  `CM_shared_peaks` = CM.shared.peaks,
                  `other_peaks` = other.peaks)

# liftover peaks from hg38 to hg19
path <- system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch <- import.chain(path)

peak.sets.hg19 <- lapply(peak.sets, function(x){unlist(liftOver(x, ch))})

# save to bed format
for(i in 1:length(peak.sets.hg19)){
  seqlevelsStyle(peak.sets.hg19[[i]]) <- "NCBI"
}

system('mkdir -p finemapping_latest_noAtrium_nonDA/annotations_for_finemapping_hg19')
system('rm finemapping_latest_noAtrium_nonDA/annotations_for_finemapping_hg19/*')

lapply(names(peak.sets.hg19), function(x){
  rtracklayer::export(peak.sets.hg19[[x]],
                      format = 'bed',
                      con = paste0('finemapping_latest_noAtrium_nonDA/annotations_for_finemapping_hg19/', x, '_hg19.bed'))})

system('cp finemapping/annotations_for_finemapping_hg19/all_eqtls_hg19.bed finemapping_latest_noAtrium_nonDA/annotations_for_finemapping_hg19/')
system('cp finemapping/annotations_for_finemapping_hg19/Coding_UCSC.bed finemapping_latest_noAtrium_nonDA/annotations_for_finemapping_hg19/')
system('cp finemapping/annotations_for_finemapping_hg19/Conserved_LindbladToh.bed finemapping_latest_noAtrium_nonDA/annotations_for_finemapping_hg19/')

# Run TORUS to get SNP prior ----------------------------------------------------------------------------
annotation_bed_files <- list.files('finemapping_latest_noAtrium_nonDA/annotations_for_finemapping_hg19',
                                   pattern = '*.bed', full.names = T)

print(annotation_bed_files)

torus.files <- prepare_torus_input_files(gwas.sumstats,
                                         annotation_bed_files,
                                         torus_input_dir = 'finemapping_latest_noAtrium_nonDA/torus_input')

torus.result <- run_torus(torus.files$torus_annot_file,
                          torus.files$torus_zscore_file,
                          option = "est-prior",
                          torus_path = "torus") # set the path to 'torus' executable.
torus.enrich <- torus.result$enrich
torus.prior <- torus.result$snp_prior

torus.result$fdr <- run_torus(torus.files$torus_annot_file,
                              torus.files$torus_zscore_file,
                              option = "fdr",
                              torus_path = "torus")$fdr
torus.fdr <- torus.result$fdr

summary(torus.result)

saveRDS(torus.result, 'finemapping_latest_noAtrium_nonDA/Torus_Enrichment_Results_Joint.rds')

torus.result <- readRDS('finemapping_latest_noAtrium_nonDA/Torus_Enrichment_Results_Joint.rds')
torus.prior <- torus.result$snp_prior

system('mkdir -p finemapping_latest_noAtrium_nonDA_newSusieR')

# finemapping using SuSiE with L = 1 ----------------------------------------------------------------

# GWAS sample size
n = 1030836

# finemapping only significant GWAS Loci (P < 5e-8)
sig.loci <- gwas.sumstats %>%
  group_by(locus) %>%
  summarise(max_mlogP = max(pval)) %>%
  filter(max_mlogP > -log10(5e-8)) %>% pull(locus)

gwas.sumstats.sigloci <- gwas.sumstats[gwas.sumstats$locus %in% sig.loci, ]
cat(length(unique(gwas.sumstats.sigloci$locus)), "GWAS significant loci. \n")

sumstats_finemap <- prepare_susie_data_with_torus_result(sumstats = gwas.sumstats.sigloci,
                                                         torus_prior = torus.prior)

sumstats_finemap <- sumstats_finemap %>% dplyr::rename(old_torus_prior = torus_pip)
cat("Finemap",length(unique(sumstats_finemap$locus)), "loci.\n")
saveRDS(sumstats_finemap, 'finemapping_latest_noAtrium_nonDA/sumstats.sigloci.rds')

# # finemapping only significant GWAS Loci (FDR < 0.1)
# sumstats_finemap_fdr <- prepare_susie_data_with_torus_result(sumstats = gwas.sumstats,
#                                                              torus_prior = torus.prior,
#                                                              torus_fdr = torus.fdr,
#                                                              fdr_thresh = 0.1)
#
# sumstats_finemap_fdr <- sumstats_finemap_fdr %>% dplyr::rename(old_torus_prior = torus_pip)
# cat("Finemap",length(unique(sumstats_finemap_fdr$locus)), "loci.\n")
# saveRDS(sumstats_finemap_fdr, 'finemapping_latest_noAtrium/sumstats.torusfdr0.1.rds')

# Run SuSiE with L = 1
bigSNP <- bigsnpr::snp_attach(rdsfile = '/project2/xinhe/1kg/bigsnpr/EUR_variable_1kg.rds')

# susie_finemap_L1 is a list of SuSiE results, one for each chunk/LD block.
susie_finemap_L1 <- run_finemapping(sumstats_finemap, bigSNP, n = n, priortype = 'torus', L = 1)

# add susie PIP information to GWAS summary stats
gwas_torus_finemapped <- merge_susie_sumstats(susie_results = susie_finemap_L1, sumstats = sumstats_finemap)
saveRDS(gwas_torus_finemapped, 'finemapping_latest_noAtrium_nonDA/AF_finemapping_result_torusprior_122loci.rds')

# Session Info -------------------
sessionInfo()

