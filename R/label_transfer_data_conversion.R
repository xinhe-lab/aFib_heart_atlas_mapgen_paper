library(Seurat)
library(SeuratDisk)

Convert("/project2/xinhe/xsun/heart_atlas/2.label_transfer/data/global_raw.h5ad", dest = "global_raw_h5seurat", overwrite = TRUE)

global_raw <- LoadH5Seurat("global_raw.h5seurat")
