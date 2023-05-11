library(Seurat)
library(harmony)
library(DoubletFinder)
library(ggplot2)
library(RColorBrewer)

ggClean <- function(){
  theme_bw() + 
    theme(text = element_text(size=18),
          legend.position = "right",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=1))
}

setwd('/project2/gca/aselewa/heart_atlas_project/')

# Load raw seurat object. Filter nUMI, and remove doublets per individual (little to no batch effects within individuals)
heart_rna <- readRDS('seurat/Heart_Seurat_RNA_all_samples_raw_noAtrium.RDS')
heart_rna <- base::subset(heart_rna, nCount_RNA > 1000 & nCount_RNA < 10000 & percent.mito < 10)

heart_rna_dataset <- SplitObject(heart_rna, split.by = 'dataset')
for(i in 1:length(heart_rna_dataset)){
  heart_rna_dataset[[i]] <- NormalizeData(object = heart_rna_dataset[[i]])
  heart_rna_dataset[[i]] <- FindVariableFeatures(heart_rna_dataset[[i]], selection.method = "vst", 
                                              nfeatures = 2000, verbose = T)
  heart_rna_dataset[[i]] <- ScaleData(heart_rna_dataset[[i]])
  heart_rna_dataset[[i]] <- RunPCA(heart_rna_dataset[[i]], verbose = T, npcs = 100)
  
  heart_rna_dataset[[i]] <- doubletFinder_v3(seu = heart_rna_dataset[[i]], 
                                          PCs = 1:10, 
                                          pN = 0.015, 
                                          pK = 0.005, 
                                          nExp = round(0.15*nrow(heart_rna_dataset[[i]]@meta.data)))
  n <- colnames(heart_rna_dataset[[i]]@meta.data)
  doubletCol <- n[startsWith("DF.classifications", x = n)]
  heart_rna_dataset[[i]]$doubletStatus <- heart_rna_dataset[[i]][[doubletCol]]
  heart_rna_dataset[[i]] <- base::subset(heart_rna_dataset[[i]], doubletStatus == "Singlet")
}

first <- heart_rna_dataset[[1]]
heart_rna_dataset <- heart_rna_dataset[-1]
heart_rna_filtered <- merge(x = first, y = heart_rna_dataset)

saveRDS(heart_rna_filtered, file = 'seurat/Heart_Seurat_RNA_all_samples_Filtered_NoAtrium.rds')

heart_rna_filtered <- NormalizeData(heart_rna_filtered, verbose = T)
heart_rna_filtered <- FindVariableFeatures(heart_rna_filtered, selection.method = "vst", 
                                             nfeatures = 2000, verbose = T)
heart_rna_filtered <- ScaleData(heart_rna_filtered)
heart_rna_filtered <- RunPCA(heart_rna_filtered, verbose = T)

heart_rna_filtered <- RunHarmony(heart_rna_filtered, group.by.vars = 'individual')

heart_rna_filtered <- RunUMAP(object = heart_rna_filtered, dims = 1:30, min.dist = 0.4, reduction = 'harmony')
heart_rna_filtered <- FindNeighbors(heart_rna_filtered, dims = 1:30, reduction ='harmony')
heart_rna_filtered <- FindClusters(object = heart_rna_filtered, resolution=0.2)

DimPlot(heart_rna_filtered, group.by = 'individual', pt.size = 0.01) 
DimPlot(heart_rna_filtered, label=T) 

# drop problematic clusters
heart_rna_filtered <- base::subset(heart_rna_filtered, seurat_clusters %in% c(0:8))
# rerun umap
heart_rna_filtered <- RunUMAP(object = heart_rna_filtered, dims = 1:30, min.dist = 0.4, reduction = 'harmony')

# visualize
DimPlot(heart_rna_filtered, label=T) 

markers <- c("TNNT2","MYBPC3","MYH7","NPPA","RGS5","ABCC9","MYH11","TAGLN","DCN","PDGFRA","PECAM1","VWF","PLP1","CD8A","LCK","CD14","FOLR2")
p <- FeaturePlot(heart_rna_filtered, features = markers)
ggsave(filename = paste0('seurat/plots/umap_combined_genes.png'), plot = p, width = 20, height=12, dpi = 300)

heart_rna_filtered <- RenameIdents(heart_rna_filtered, 
                                       '0'='Endothelial', 
                                       '1'='Fibroblast',
                                       '2'='Cardiomyocyte',
                                       '3'='Pericyte',
                                       '4'='Myeloid',
                                       '5'='Smooth Muscle',
                                       '6'='Endothelial',
                                       '7'='Lymphoid',
                                       '8'='Neuronal')

palette <- readRDS('notebooks/palette.rds')

p <- DimPlot(heart_rna_filtered, pt.size = 0.05, label = T, cols=palette) + ggtitle("scRNA-seq - 67114 Nuclei") + theme_set(theme_gray()) + ggClean() 
ggsave(filename = paste0('seurat/plots/umap_combined_celltypes.png'), plot = p, width = 8, height=6, dpi = 300)

saveRDS(heart_rna_filtered, file = "seurat/Heart_RNA_Processed_Combined_NoAtrium.rds")

markers <- FindAllMarkers(object = srna, logfc.threshold = 0.25, only.pos = T)
saveRDS(markers, file = 'seurat/CM_rna_markers.rds')

# number of cells per region by donor
stats.df <- data.frame(numi=srna$nCount_RNA, donor=srna$individual, region=srna$region) %>% 
    group_by(region, donor) %>% 
    summarise(ncell=n())
pdf('manuscript_figures/nNuc_RNA_postQC.pdf',width=10,height=8)
ggplot(stats.df, aes(x = region, y = ncell, fill=donor)) + geom_bar(stat="identity", position="dodge") + ggClean(rotate_axis = T) + ylab("Number of Nuclei") + xlab("") + scale_fill_brewer(palette = "Set2")
dev.off()

pdf('manuscript_figures/RNA_nUMI_UMAP.pdf', width=10, height=8)
FeaturePlot(srna, features = 'nCount_RNA') + scale_color_gradientn(colours = c("royalblue1","white","red")) + xlab('UMAP1') +  ylab('UMAP2')
dev.off()

pdf('manuscript_figures/RNA_Genes_UMAP.pdf', width=10, height=8)
FeaturePlot(srna, features = 'nFeature_RNA') + scale_color_gradientn(colours = c("royalblue1","white","red")) + xlab('UMAP1') + ylab('UMAP2')
dev.off()


pdf('manuscript_figures/RNA_celltypes_umap.pdf',width=8, height=6)
custom_dim_plot(seurat = srna, palette = palette, label = T, legend=F) + ggtitle('snATAC-seq - Transferred Label')
dev.off()
