library(finemappeR)
library(liftOver)
library(ComplexHeatmap)
source('R/analysis_utils.R')
setwd('/project2/gca/aselewa/heart_atlas_project/')
celltype_ideal_order <- c("Cardiomyocyte","Smooth Muscle","Pericyte","Endothelial","Fibroblast","Neuronal", "Lymphoid","Myeloid")

gwas.sumstats <- readRDS('GWAS/summary_statistics/aFib/ebi-a-GCST006414_aFib.df.rds')

satac <- loadArchRProject('ArchR/ArchR_heart_latest_noAtrium/')
motif.pos <- getPositions(satac)
union.set <- getPeakSet(satac)

tfs.summary <- read_tsv('TF_Correlation_Results/Final_TFs_summary_06272021.tsv')

motifPositions <- getPositions(satac)
#markerMotifs <- unlist(lapply(c("ESRRB","ESRRG","MEF2A","MEF2D","TBX5","TBX2","GATA4","GATA5"), function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs <- unlist(lapply(c("TBX5","ETV1","PITX2","PRRX1","NKX25","GATA4"), function(x) grep(x, names(motifPositions), value = TRUE)))
motifPositions <- motifPositions[markerMotifs]

#calculate footprints for each TF
seFoot <- getFootprints(
  ArchRProj = satac, 
  positions = motifPositions, 
  groupBy = "CellTypes"
)

# plot footprints
plotFootprints(
  seFoot = seFoot,
  ArchRProj = satac, 
  normMethod = "Divide",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)


motifPositions <- motifPositions[-3]

motifPositions.in.OCRs <- lapply(motifPositions, function(x){subsetByOverlaps(x, union.set, minoverlap = 2)})
motifPositions.in.OCRs.hg19 <- lapply(motifPositions.in.OCRs, function(x){hg38ToHg19(x)})

for(i in 1:length(motifPositions.in.OCRs.hg19)){
  seqlevelsStyle(motifPositions.in.OCRs.hg19[[i]]) <- "NCBI"
  motifPositions.in.OCRs.hg19[[i]] <- motifPositions.in.OCRs.hg19[[i]] + 50
}

merged.OCRs.hg19 <- unlist(GRangesList(motifPositions.in.OCRs.hg19))

system('mkdir GWAS/footprinting_annotations')
rtracklayer::export(merged.OCRs.hg19, con = 'GWAS/footprinting_annotations/Merged_Motifs_50bp.bed', format = 'bed')

# lapply(names(motifPositions.in.OCRs.hg19), function(x){
#   rtracklayer::export(motifPositions.in.OCRs.hg19[[x]], con = paste0('GWAS/footprinting_annotations/',x,'.bed'), format = 'bed')
# })

union.set.hg19 <- hg38ToHg19(union.set)
seqlevelsStyle(union.set.hg19) <- "NCBI"
rtracklayer::export(union.set.hg19, con = 'GWAS/footprinting_annotations/union_OCRs.bed', format = 'bed')


annotations <- list.files(path = 'GWAS/footprinting_annotations', pattern = '.bed', full.names = T)
enrich.res <- list()
for(w in seq_along(annotations)){
  
  cleaned.gwas.annots <- annotator(gwas.sumstats, annotations = annotations[w])
  
  readr::write_tsv(x = cleaned.gwas.annots[,-c(1:6,8:12)], path = 'torus_annotations.txt.gz', col_names = T)
  readr::write_tsv(x = gwas.sumstats[,c('snp','locus','zscore')], path = 'torus_zscore.txt.gz', col_names = T)
  
  torus.result <- RunTorus(torus_annot_file = 'torus_annotations.txt.gz', torus_zscore_file = 'torus_zscore.txt.gz') 
  
  enrich.res[[basename(annotations[w])]] <- torus.result$enrich
}


saveRDS(enrich.res, 'GWAS/Torus_TF_footprint_enrichment_50bp_extend.rds')

enrich.df <- Reduce(rbind, enrich.res)
enrich.df <- enrich.df[enrich.df$term != "Intercept",]
enrich.df$term <- sub('_.*','', enrich.df$term)
enrich.df

ggplot(enrich.df, aes(y = estimate/log(2), x = term)) + 
  geom_point() + 
  geom_errorbar(mapping = aes(ymin = low/log(2), ymax=high/log(2))) + 
  coord_flip() +
  ggClean() +
  geom_hline(yintercept = 1, linetype='dashed') +
  xlab('TF') + ylab('Log2 Fold-of-Enrichment') +
  ggtitle('Footprint Enrichment in AFib GWAS (+/- 50bp)')

