library(ArchR)
library(ggplot2)
library(tidyverse)

addArchRThreads(threads = 12) 
setwd('/project2/gca/aselewa/heart_atlas_project/')
macs2 <- '/project2/gca/software/miniconda3/bin/macs2'

source('R/analysis_utils.R')

archr_project_path <- 'ArchR/ArchR_heart_latest_noAtrium/'

projHeart <- loadArchRProject(archr_project_path)

projHeart <- addGroupCoverages(ArchRProj = projHeart, groupBy = "CellTypes", force = T, maxCells = 10000)
projHeart <- addReproduciblePeakSet(ArchRProj = projHeart, groupBy = "CellTypes", pathToMacs2 = macs2, cutOff = 0.01, verbose = T)
projHeart <- addPeakMatrix(projHeart, force = T)

# cell-type specific peaks
markersPeaks <- getMarkerFeatures(ArchRProj = projHeart, 
                                  useMatrix = "PeakMatrix", 
                                  groupBy = "CellTypes", 
                                  bias = c("TSSEnrichment", "log10(nFrags)"))

saveRDS(markersPeaks, paste0(archr_project_path,'/PeakCalls/DA_markerPeaks.rds'))

markers <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.1", returnGR = T)
saveRDS(markers, file = paste0(archr_project_path,'/PeakCalls/DA_MARKERS_FDRP_10_log2FC_0.rds'))

saveArchRProject(projHeart) 

# Motif Enrichment
projHeart <- addMotifAnnotations(projHeart, name = "Motif")
projHeart <- addDeviationsMatrix(ArchRProj = projHeart, peakAnnotation = "Motif", force = T)

# Co-accessibility
satac <- addCoAccessibility(ArchRProj = satac, reducedDims = 'harmony', maxDist = 1e6)

# BigWigs by cell-type
getGroupBW(ArchRProj = satac, groupBy = "CellTypes")

saveArchRProject(satac) 




