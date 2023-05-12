library(Seurat)
library(Azimuth)  ###module load hdf5_hl
library(SeuratData)
library(patchwork)


dat_new <- readRDS("/project2/xinhe/xsun/heart_atlas/2.label_transfer/data/global_raw.rds")
dat_our <- readRDS("/project2/gca/aselewa/heart_atlas_project/seurat/Heart_RNA_Processed_Combined_NoAtrium.rds")
dat_our$celltypes <- Idents(dat_our)

dat_ref <- AzimuthReference(pbmc_small,refUMAP = "umap",
                            refDR = "spca",
                            refAssay = "RNA",
                            dims = 1:50,
                            k.param = 31,
                            plotref = "umap",
                            plot.metadata = NULL,
                            ori.index = NULL,
                            colormap = NULL,
                            assays = NULL,
                            metadata = NULL,
                            reference.version = "0.0.0",
                            verbose = FALSE)

