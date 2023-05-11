
library(ggplot2)
library(dplyr)
source('R/analysis_utils.R')

rna.meta <- readRDS('seurat/RNA_seq_metadata.rds')
atac.meta <- readRDS('ATAC_metadata_df.rds')

atac.meta$sample <- rna.meta$sample
comb.meta <- rbind(rna.meta, atac.meta)
comb.meta$Assay <- c(rep("RNA",18),rep("ATAC",18))

collapsed.comb.meta <- comb.meta %>% 
  group_by(individual, Assay) %>% 
  summarise(num_cells=sum(ncells), mean_saturation=mean(saturation, na.rm=T), sd_saturation=sd(saturation, na.rm=T))

p1 <- ggplot(collapsed.comb.meta, aes(x=individual,y=mean_saturation, fill=Assay)) + 
  geom_bar(position="dodge", stat="identity") + 
  ggClean(rotate_axis=TRUE) + xlab("Donor") + ylab("% Saturation") + 
  geom_errorbar(aes(ymin=mean_saturation-sd_saturation/2, ymax=mean_saturation+sd_saturation/2), width=.2,
                position=position_dodge(.9))

p2 <- ggplot(collapsed.comb.meta, aes(x=individual,y=num_cells, fill=Assay)) + 
  geom_bar(position="dodge", stat="identity") + 
  ggClean(rotate_axis=TRUE) + xlab("Donor") + ylab("Cells Detected") 

ggsave(filename = "individual_saturation_qc.png", p1, dpi=300, width=8, height=6)
ggsave(filename = "individual_cells_detected.png", p2, dpi=300, width=8, height=6)

individuals <- unique(comb.meta$individual)
for(ind in individuals){
  curr_df <- comb.meta[comb.meta$individual==ind,]
  p1 <- ggplot(curr_df, aes(x=regions,y=saturation, fill=Assay)) + 
    geom_bar(position="dodge", stat="identity") + 
    ggClean(rotate_axis=TRUE) + xlab("Region") + ylab("% Saturation") + ggtitle(ind) + ylim(c(0, 100))
  
  p2 <- ggplot(curr_df, aes(x=regions,y=ncells, fill=Assay)) + 
    geom_bar(position="dodge", stat="identity") + 
    ggClean(rotate_axis=TRUE) + xlab("Region") + ylab("Cells Detected") + ggtitle(ind) + ylim(c(0,8000))
  
  ggsave(filename = paste0("regions_",ind,"_saturation_qc.png"), p1, dpi=300, width=8, height=6)
  ggsave(filename = paste0("regions_",ind,"_ncells_qc.png"), p2, dpi=300, width=8, height=6)
}


