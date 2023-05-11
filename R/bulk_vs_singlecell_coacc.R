library(GenomicRanges)
library(ggplot2)
setwd('/project2/gca/aselewa/heart_atlas_project/')
source('R/analysis_utils.R')

coacc.list <- readRDS('ArchR/ArchR_heart_latest_noAtrium/CoAccessibility/Coaccess_celltype_list.rds')
CM.coacc <- coacc.list$Cardiomyocyte
CM.coacc <- CM.coacc[CM.coacc$correlation > 0.4,]

sc.pairs <- bind_rows(coacc.list)
sc.pairs <- sc.pairs[sc.pairs$correlation>0.2,]

bulk.coacc <- readRDS('ArchR/ArchR_heart_latest_noAtrium/CoAccessibility/bulk_coacc_corr0.df.rds')
bulk.coacc <- bulk.coacc[bulk.coacc$correlation > 0.5,]

CM.pairs <- paste0(CM.coacc$queryHits,'_', CM.coacc$subjectHits)
bulk.pairs <- paste0(bulk.coacc$queryHits,'_', bulk.coacc$subjectHits)
sc.pairs.str <- paste0(sc.pairs$queryHits,'_',sc.pairs$subjectHits)

length(intersect(CM.pairs, bulk.pairs))/length(CM.pairs)

length(intersect(CM.pairs, bulk.pairs))/length(bulk.pairs)

length(intersect(sc.pairs.str, bulk.pairs))/length(bulk.pairs)

