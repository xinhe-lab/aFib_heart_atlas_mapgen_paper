setwd('/project2/gca/aselewa/heart_atlas_project/')

library(ArchR)
source('R/analysis_utils.R')

ATAC_SAMPLES <- c("MW200804AA","Deep_MW200804AC","Deep_MW200804AD","SP-HE-MW200928E1ATAC-413ATAC",
                  "SP-HE-HE200915ATAC-175ATAC","Deep_SP-HE200915ATAC-360ATAC","Deep_SP-HE200915ATAC-397ATAC","SP-HE-MW200928E1ATAC-398ATAC",
                  "SP-HE-MW200928E2ATAC-175ATAC","Deep_SP-MW200928E2ATAC-367ATAC","Deep_SP-MW200928E2ATAC-407ATAC","SP-HE-MW200928E1ATAC-408ATAC")
ATAC_INDIVIDUALS <- c(rep("02207",4),rep("02336",4),rep("03231",4))
ATAC_REGIONS <- rep(c("Septum","Right Ventricle","Left Ventricle","Apex"), 3)

palette <- readRDS('notebooks/palette.rds')

project_path <- 'ArchR/ArchR_heart_latest_noAtrium_latest/'

addArchRThreads(threads = 5)
addArchRGenome("hg38")
archr_project <- ArchRProject(
  ArrowFiles = list.files(path = 'ArchR/ArchR_ArrowFiles_5kmin_6TSSmin/', pattern = '*.arrow', full.names = T),
  outputDirectory = project_path,
  copyArrows = T
)
archr_project$regions <- plyr::mapvalues(x = archr_project$Sample, from = ATAC_SAMPLES, to = ATAC_REGIONS)
archr_project$individual <- plyr::mapvalues(x = archr_project$Sample, from = ATAC_SAMPLES, to = ATAC_INDIVIDUALS)
saveArchRProject(archr_project)

archr_project <- filterDoublets(ArchRProj = archr_project)

archr_project <- addIterativeLSI(ArchRProj = archr_project, 
                                 useMatrix = "TileMatrix", 
                                 name = "IterativeLSI", 
                                 iterations = 2, 
                                 varFeatures = 20000, 
                                 force = TRUE)

archr_project <- addHarmony(ArchRProj = archr_project, 
                            reducedDims = "IterativeLSI", 
                            name = "harmony", 
                            groupBy = "individual", 
                            force = TRUE)

archr_project <- addClusters(input = archr_project, 
                             reducedDims =  "harmony", 
                             resolution = 0.8, 
                             force = TRUE)

archr_project <- addUMAP(ArchRProj = archr_project, 
                         reducedDims = "harmony", 
                         force=TRUE, 
                         minDist = 0.6)

archr_project <- saveArchRProject(archr_project, load = T)

# check clustering-visualization agreement
custom_archr_plot(archr_project = archr_project, group.by="Clusters", pt.size=0.3, alpha=0.7, label=T, legend = T)

# remove clusters that aren't well defined (doublets or low quality cells that passed QC)
archr_project <- subsetCells(ArchRProj = archr_project, 
                             cellNames = archr_project$cellNames[!(archr_project$Clusters %in% c("C2","C13","C14","C15"))])

archr_project <- addUMAP(ArchRProj = archr_project, 
                         reducedDims = "harmony", 
                         force=TRUE, 
                         minDist = 0.6)

archr_project <- addImputeWeights(archr_project)
markers <- c("TNNT2","MYBPC3","MYH7","NPPA","RGS5","ABCC9","MYH11","TAGLN","DCN","PDGFRA","PECAM1","VWF","PLP1","CD8A","LCK","CD14","FOLR2")
p <- plotEmbedding(ArchRProj = archr_project, colorBy = "GeneScoreMatrix", name = markers, embedding = "UMAP")
plotPDF(plotList = p, name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", ArchRProj = archr_project, addDOC = FALSE, width = 5, height = 5)

p <- custom_archr_plot(archr_project = archr_project, pt.size=0.3, alpha=0.7, label=F, legend = T)
ggsave(filename = paste0(project_path,"/Plots/umap-regions.png"), plot = p, dpi=150, width=8, height=6)

p <- custom_archr_plot(archr_project = archr_project, pt.size=0.3, alpha=0.7, label=F, legend = T, group.by = 'individual')
ggsave(filename = paste0(project_path,"/Plots/umap-donors.png"), plot = p, dpi=150, width=8, height=6)

p <- custom_archr_plot(archr_project = archr_project, group.by="Clusters", pt.size=0.3, alpha=0.7, label=T, legend = T)
ggsave(filename = paste0(project_path,"/Plots/umap-Clusters.png"), plot = p, dpi=150, width=8, height=6)

p <- getQCPlots(archr_project)
ggsave(filename = paste0(project_path,"/Plots/VlnQCPlots.png"), plot = p, dpi=150, width=14, height=10)

# map cluster IDs to cell-types
cluster.ids <- paste0("C",1:19)
new.ids <- c(rep("Cardiomyocyte", 6),"Neuronal","Smooth Muscle","Pericyte",rep("Fibroblast",3),rep("Doublet",3),rep("Endothelial",2),"Myeloid","Lymphoid")
archr_project$CellTypes <- RenameIdentity(idents = archr_project$Clusters, from = cluster.ids, to = new.ids)

p <- custom_archr_plot(archr_project = archr_project, group.by="CellTypes", palette = palette, pt.size=0.3, alpha=0.7, label=T, legend = T)
ggsave(filename = "ArchR/ArchR_heart_latest_noAtrium/Plots/umap-CellTypes.png", plot = p, dpi=150, width=8, height=6)
saveArchRProject(archr_project)

# better heatmap of gene activities
markers <- c("TNNT2","RGS5","TAGLN","DCN","VWF","PLP1","CD8A","CD14")
p <- plotEmbedding(ArchRProj = archr_project, colorBy = "GeneScoreMatrix", name = markers,embedding = "UMAP")

plist <- list()
z <- names(p)
for(i in 1:length(z)){
  plist[[i]] <- p[[i]] + ggClean() + 
    LegendOff() + 
    xlab("") + 
    ylab("") + 
    ggtitle(z[i]) + 
    theme(text = element_text(size=12)) + scale_fill_gradientn(colours = c("lightblue","yellow","red"))
}
do.call("grid.arrange", plist)

# number of nuclei in the QC version
stats.df <- data.frame(nfrag=archr_project$nFrags, donor=archr_project$individual, region=archr_project$regions) %>% 
    group_by(region, donor) %>% 
    summarise(ncell=n())
pdf('manuscript_figures/nNuc_postQC.pdf',width=10,height=8)
ggplot(stats.df, aes(x = region, y = ncell, fill=donor)) + geom_bar(stat="identity", position="dodge") + ggClean(rotate_axis = T) + ylab("Number of Nuclei") + xlab("") + scale_fill_brewer(palette = "Set2")
dev.off()

# visualize regions
pdf('manuscript_figures/ATAC_regions_umap.pdf',width=10, height=8)
custom_archr_plot(archr_project = satac, group.by = 'regions', palette = c("purple","green","blue","red"), pt.size = 0.7, alpha = 0.5)
dev.off()

pdf('manuscript_figures/ATAC_nfrags.pdf',width=10, height=8)
plotEmbedding(satac, colorBy = "cellColData", name = "nFrags") + ggClean() + xlab('UMAP1') + ylab('UMAP2')
dev.off()

pdf('manuscript_figures/ATAC_doublet_enrichment.pdf',width=10, height=8)
plotEmbedding(satac, colorBy = "cellColData", name = "DoubletEnrichment") + ggClean() + xlab('UMAP1') + ylab('UMAP2')
dev.off()

#integration with RNA
srna <- readRDS('seurat/Heart_RNA_Processed_Combined_NoAtrium.rds')
srna$celltypes <- Idents(srna)
archr_project <- addGeneIntegrationMatrix(
    ArchRProj = archr_project,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "harmony",
    seRNA = srna,
    addToArrow = FALSE,
    groupRNA = "celltypes",
    nameCell = "predictedCell_srna",
    nameGroup = "predictedGroup_srna",
    nameScore = "predictedScore_srna",
    force =TRUE
)
saveArchRProject(archr_project)

pdf('manuscript_figures/ATAC_transferred_label_umap.pdf',width=8, height=6)
custom_archr_plot(archr_project = satac, 
                  group.by = 'predictedGroup_srna', 
                  palette = palette, 
                  legend = F,
                  label = T) + ggtitle('snATAC-seq - Transferred Label')
dev.off()

pdf('manuscript_figures/ATAC_umap.pdf',width=8, height=6)
custom_archr_plot(archr_project = satac, 
                  group.by = 'CellTypes', 
                  palette = palette, 
                  legend = F,
                  label = T) + ggtitle('snATAC-seq - Transferred Label')
dev.off()
