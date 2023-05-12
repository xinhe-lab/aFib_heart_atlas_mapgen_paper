library(SCAVENGE)
library(gchromVAR)
library(chromVAR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SummarizedExperiment)
library(data.table)
library(dplyr)



### peak matrix
#projHeart <- loadArchRProject("/project2/xinhe/xsun/heart_atlas/3.ATAC_peak/CopyOfArchR_Heart_Latest_Backup/")
#peakmatrix <- getMatrixFromProject(projHeart,useMatrix = "PeakMatrix")
#saveRDS(peakmatrix, file = "peakmatrix.rds")
peakmatrix <- readRDS("peakmatrix.rds")
peakmatrix <- addGCBias(peakmatrix, genome = BSgenome.Hsapiens.UCSC.hg38)
assayNames(peakmatrix) <- "counts"
peakmatrix_bg <- getBackgroundPeaks(peakmatrix, niterations=200)

### trait file
trait_file <- "finemapping_aFib_hg38_new.bed"   #####
trait_import <- importBedScore(rowRanges(peakmatrix), trait_file, colidx=5)
peakmatrix_DEV <- computeWeightedDeviations(peakmatrix, trait_import,background_peaks = peakmatrix_bg)

z_score_mat <- data.frame(colData(peakmatrix), z_score=t(assays(peakmatrix_DEV)[["z"]]) %>% c)

#####
seed_idx <- seedindex(z_score_mat$z_score, 0.05)

# Cells with enriched P < 0.05: 2382
# Percent: 8.92%
# The top 5% of cells (N=1335) were selected as seed cells

scale_factor <- cal_scalefactor(z_score=z_score_mat$z_score, 0.01)
# Scale factor is calculating from most enriched 1% of cells
##

peak_by_cell_mat <- assay(peakmatrix)
tfidf_mat <- tfidf(bmat=peak_by_cell_mat, mat_binary=TRUE, TF=TRUE, log_TF=TRUE)

####
lsi_mat <- do_lsi(tfidf_mat, dims=30)
mutualknn30 <- getmutualknn(lsi_mat, 30)

######
np_score <- randomWalk_sparse(intM=mutualknn30, rownames(mutualknn30)[seed_idx], gamma=0.05)

# Stationary step: 151
# Stationary Delta: 9.91057731496478e-06

omit_idx <- np_score==0
sum(omit_idx)
# 492

mutualknn30 <- mutualknn30[!omit_idx, !omit_idx]
np_score <- np_score[!omit_idx]
TRS <- np_score %>% capOutlierQuantile(., 0.95) %>% max_min_scale
TRS <- TRS * scale_factor
cellcoldata <- data.frame(z_score_mat[!omit_idx, ], seed_idx[!omit_idx], np_score, TRS)

save(cellcoldata, file = "TRS_result_new.rdata")




######plots





