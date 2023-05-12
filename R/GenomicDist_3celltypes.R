library(GenomicDistributions)
library(GenomicDistributionsData)  #### run in R/4.2.0
library(ggplot2)
library(gridExtra)
library(rtracklayer)

setwd("/project2/xinhe/xsun/heart_atlas/1.OCR_compare")

celltype_our <- c("Cardiomyocyte","Endothelial","Fibroblast")

celltype_new <- c("Ventricular_cardiomyocyte","Endothelial","Fibroblast")
celltype_new_name <- c("vCardiomyocyte","Endothelial","Fibroblast")


p<-list()
for (i in 1:length(celltype_our)) {
  
  file_our <- paste0("/project2/gca/aselewa/heart_atlas_project/ArchR/ArchR_heart_latest_noAtrium/PeakCalls/",celltype_our[i],"-reproduciblePeaks.gr.rds")
  dat_our <- readRDS(file_our)
  
  extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                            qValue = "numeric", peak = "integer")
  dat_new <- import(paste0("/project2/xinhe/xsun/heart_atlas/1.OCR_compare/peakdata_ren/celltype/",celltype_new[i],".narrowPeak.bed"), format = "BED",
                    extraCols = extraCols_narrowPeak)
  dat_new <- unique(dat_new)
  queryList <- GRangesList("Hocker's"=dat_new, "Ours"=dat_our)
  cal <- calcPartitionsRef(queryList, "hg38")
  
  p[[i]] <- plotPartitions(cal) +
    ggtitle(paste0(celltype_our[i],"(this study)","\n and \n ",celltype_new_name[i], "(Hocker et al.)")) +
    theme(plot.title = element_text(size=35),
          legend.text=element_text(size=30),
          legend.title=element_text(size=30)) +
    theme(axis.title.x = element_text(size = 35),
          axis.text.x = element_text(size = 30, color = "black"),
          axis.title.y = element_text(size = 35),
          axis.text.y = element_text(size = 30, color = "black")) 
  
}

all <- grid.arrange(p[[1]],p[[2]],p[[3]], nrow = 1)
#ggsave(filename = "/project2/xinhe/xsun/atlas/genomicdist_3celltypes_larger.jpg", plot = all, dpi = 600,width = 30, height = 10, limitsize = FALSE )
ggsave(filename = "plots/genomicdist_3celltypes_larger.pdf", plot = all, device = "pdf",dpi = 600,width = 30, height = 10, limitsize = FALSE )

