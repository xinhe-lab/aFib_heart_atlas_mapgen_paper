library(genomation)
library(GenomicRanges)
library(ComplexHeatmap)
library(ggplot2)


setwd("/project2/xinhe/xsun/heart_atlas/5.OR/")

load("/project2/xinhe/xsun/heart_atlas/5.OR/data/OCR_celltypespecific.rdata")
load("/project2/xinhe/xsun/heart_atlas/5.OR/data/encode.beds.rdata")


subsetByOverlapProp <- function(q, s, minPoverlap, maxgap=0){
  
  hits <- GenomicRanges::findOverlaps(query = q, subject = s, maxgap = maxgap)
  overlaps <- pintersect(q[queryHits(hits)], s[subjectHits(hits)])
  percentOverlap <- width(overlaps) / width(q[queryHits(hits)])
  hits <- hits[percentOverlap >= minPoverlap]
  
  return(hits)
}


odd_all <- c()
for (i in 1:length(encode.beds)) {
  
  gr_bulk <- encode.beds[[i]]
  
  odd1 <- c()
  odd2 <- c()
  for (j in 1:length(markers)) {
    gr_celltype <- markers[[j]]
    
    hit_overlap_with_bulk <- subsetByOverlapProp(q = gr_celltype, s = gr_bulk, minPoverlap = 0.1)
    overlap_with_bulk <- gr_celltype[unique(queryHits(hit_overlap_with_bulk)),]
    
    ### count length of overlap for each random data set
    overlap_with_random_tmp <- c()
    for (t in 1:10) {
      gr_random <- readBed(paste0("/project2/xinhe/xsun/heart_atlas/5.OR/data/data_shuffled_samesize_10times_exclude_itself/",names(encode.beds)[i],"_",t,"_shuffled.bed" ))
      
      hit_overlap_with_random <- subsetByOverlapProp(q = gr_celltype, s = gr_random, minPoverlap = 0.1)
      overlap_with_random_tmp[t] <- length(gr_celltype[unique(queryHits(hit_overlap_with_random)),])
      
    }
    
    overlap_with_random <- mean(overlap_with_random_tmp) ###length
    odd1[j] <- (length(overlap_with_bulk) / length(gr_bulk)) / ( 1- (length(overlap_with_bulk) / length(gr_bulk)))
    odd2[j] <- (overlap_with_random  / length(gr_random))  / (1 - (overlap_with_random  / length(gr_random))) 
    
  }
  odd <- odd1/odd2
  odd_all <- rbind(odd_all,odd)
  
}


colnames(odd_all) <- names(markers)
rownames(odd_all) <- names(encode.beds)

odd_re <- odd_all[,c("Cardiomyocyte","Pericyte","Endothelial","Fibroblast","Lymphoid","Myeloid","Neuronal")]
#odd_re <- odd_all[,c("Cardiomyocyte","Pericyte","Endothelial","Fibroblast","Lymphoid","Myeloid","Neuronal","Smooth Muscle")]
#odd_re <- odd_all[,c("Cardiomyocyte","Pericyte","Endothelial","Fibroblast","Lymphoid","Myeloid")]
odd_heatmap <- odd_re[c("Heart_RVLV","DermisEndothelialBlood","Tcell","Bcell","CD14Monocyte"),]

rownames(odd_heatmap) <- c("Heart \nLV+RV","Dermis \nEndothelial","T cell","B cell","CD14+ \nMonocyte")
colnames(odd_heatmap) <- c("CM","Peri","Endo","Fibro","Lymph", "Mye","Neuro")
#colnames(odd_heatmap) <- c("CM","Peri","Endo","Fibro","Lymph", "Mye")

pdf('heatmap_it_samesize_10times.pdf', width=12, height=10)
Heatmap(matrix = odd_heatmap, cluster_rows = F, 
        cluster_columns = F, show_row_names = T, row_names_side = "left",
        name = "Odd Ratio",
        col = circlize::colorRamp2(c(0,8,60), c("lightblue","white","firebrick")), 
        row_title = "ENCODE DNase",
        column_title = "Cell Types",
        row_gap = unit(1, "mm"),
        column_gap = unit(1, "mm"),
        na_col = "white",
        rect_gp = gpar(col = "black", lwd = 0.5),
        use_raster = F)
dev.off()

odd_log <- log(odd_heatmap)
pdf('heatmap_log_it_samesize_10times.pdf', width=12, height=10)
Heatmap(matrix = odd_log, cluster_rows = F, 
        cluster_columns = F, show_row_names = T, row_names_side = "left",
        name = "log(OR)",
        col = circlize::colorRamp2(c(0,2,4.1), c("lightblue","white","firebrick")), 
        row_title = "ENCODE DNase",
        column_title = "Cell Types",
        row_gap = unit(1, "mm"),
        column_gap = unit(1, "mm"),
        na_col = "white",
        rect_gp = gpar(col = "black", lwd = 0.5),
        use_raster = F)
dev.off()




#####H3K27ac

ggClean <- function(rotate_axis=FALSE){
  tm <- theme_bw() + 
    theme(text = element_text(size=18),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=1))
  if(rotate_axis){
    tm <- tm + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
  
  return(tm)
  
}

odd_re <- odd_all[,c("Cardiomyocyte","Pericyte","Endothelial","Fibroblast","Lymphoid","Myeloid","Neuronal")]
odd_h3k <- data.frame(odd_re["Heart_LVRV",])
odd_h3k <- cbind(rownames(odd_h3k),odd_h3k)
colnames(odd_h3k) <- c("celltypes","Odd_ratio")
odd_h3k$celltypes <- c("CM","Peri","Endo","Fibro","Lymph", "Mye","Neuro")
odd_h3k$celltypes <- factor(odd_h3k$celltypes,levels = odd_h3k$celltypes)

pdf('h3k27ac_hist_samesize_10times.pdf', width=4, height=4)
ggplot(odd_h3k, aes(x=celltypes, y=Odd_ratio)) + 
  geom_bar(stat='identity', width=0.7) + 
  ggClean(rotate_axis = T) + 
  ylab('Odd Ratio') + 
  xlab('')
dev.off()

odd_h3k$Odd_ratio <- log(odd_h3k$Odd_ratio)
pdf('h3k27ac_hist_log_it_samesize_10times.pdf', width=4, height=4)
ggplot(odd_h3k, aes(x=celltypes, y=Odd_ratio)) + 
  geom_bar(stat='identity', width=0.7) + 
  ggClean(rotate_axis = T) + 
  ylab('log(OR)') + 
  xlab('')
dev.off()
  
  
  
