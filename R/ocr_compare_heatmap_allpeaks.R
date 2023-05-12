library(rtracklayer)
library(ggplot2)
library(ComplexHeatmap)

setwd("/project2/xinhe/xsun/heart_atlas/1.OCR_compare")

genomicCorr.jaccard = function(query, reference, restrict = NULL) {
  if(is.null(restrict)) {
    res = sum(width(intersect(query, reference))) / sum(as.numeric(width(union(query, reference))))
  } else {
    gr1 = intersect(query, reference)
    gr1 = intersect(gr1, restrict)
    
    gr2 = union(query, reference)
    gr2 = intersect(gr2, restrict)
    res = sum(width(gr1)) / sum(width(gr2))
  }
  return(res)
}

celltype_our <- c("Cardiomyocyte","Endothelial","Fibroblast","Lymphoid","Myeloid","Neuronal","Smooth.Muscle","Pericyte")
celltype_name <- c("Cardiomyocyte","Endothelial","Fibroblast","Lymphoid","Myeloid","Neuronal","Smooth Muscle","Pericyte")

celltype_new <- c("Ventricular_cardiomyocyte","Atrial_cardiomyocyte","Endothelial","Fibroblast","Lymphocyte","Macrophage","Nervous","Smooth_muscle","Adipocyte")
celltype_new_name <-c("vCardiomyocyte","aCardiomyocyte","Endothelial","Fibroblast","Lymphocyte","Macrophage","Nervous","Smooth muscle","Adipocyte")

folder_new_data <- "/project2/xinhe/xsun/heart_atlas/1.OCR_compare/peakdata_ren/celltype/"

jcd_matrix <- matrix(,nrow=length(celltype_our), ncol = length(celltype_new))
for(i in 1:length(celltype_our)) {
  
  file_our <- paste0("/project2/gca/aselewa/heart_atlas_project/ArchR/ArchR_heart_latest_noAtrium/PeakCalls/",celltype_our[i],"-reproduciblePeaks.gr.rds")
  dat_our <- readRDS(file_our)
  for (j in 1:length(celltype_new) ) {
    
    extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                              qValue = "numeric", peak = "integer")
    files <- paste0(folder_new_data,celltype_new[j],".narrowPeak.bed")
    dat_new <- import(files, format = "BED",extraCols = extraCols_narrowPeak)
    dat_new <- unique(dat_new) ##########
    
    jcd_tmp <- genomicCorr.jaccard(dat_new,dat_our)
    
    jcd_matrix[i,j] <- jcd_tmp
    
  }
}
colnames(jcd_matrix) <- celltype_new_name
rownames(jcd_matrix) <- celltype_name

save(jcd_matrix, file = "jaccard_heatmap/jcd_matrix_all_final.rdata")

pdf('plots/jcd_matrix_all_heatmap_final.pdf', width=10, height=10)
#jpeg('/project2/xinhe/xsun/atlas/jcd_matrix_all_heatmap_uniq_final.jpg', width=1100, height=1000, res=200)
Heatmap(matrix = jcd_matrix, cluster_rows = F, 
        cluster_columns = F, show_row_names = T, row_names_side = "left",
        name = "Jaccard Index",
        col = circlize::colorRamp2(c(0, 0.45), c("white","firebrick")), ###modify the figures
        row_title = NULL,
        column_title = NULL,
        row_gap = unit(1, "mm"),
        column_gap = unit(1, "mm"),
        na_col = "white",
        use_raster = T,
        column_names_gp = grid::gpar(fontsize = 18),
        row_names_gp = grid::gpar(fontsize = 18))
dev.off()


