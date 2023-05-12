library(ArchR)
library(ggplot2)
library(parallel)


projHeart <- loadArchRProject("/project2/xinhe/xsun/heart_atlas/3.ATAC_peak/ArchR_Heart_Latest_Backup/")
load("/project2/xinhe/xsun/heart_atlas/4.scavenge/TRS_result_new.rdata")


p <- ggplot(cellcoldata, aes(x=as.factor(CellTypes), y=TRS)) +
  geom_boxplot() +
  theme_bw(base_line_size =0.3) +
  labs(x = "Celltypes", y= "SCAVENGE TRS") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14, color = "black",angle = 90, vjust = 0.5, hjust=1),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14, color = "black"))
ggsave(p, filename = "boxplot_unif_new.pdf", dpi = 600, width= 5, height = 5)



####
projHeart2 <-addCellColData( ArchRProj = projHeart, data = cellcoldata[,29], cells = rownames(cellcoldata),name = "TRS")
projHeart2$TRS[is.na(projHeart2$TRS)] <- 0

p2 <- plotEmbedding(ArchRProj = projHeart2, colorBy = "cellColData", name = "TRS", embedding = "UMAP")
#ggsave(p2, filename = "UMP_TRS.jpg", dpi = 600, width= 5, height = 5)
ggsave(p2, filename = "UMP_TRS_new.pdf", dpi = 600, width= 5, height = 5)
