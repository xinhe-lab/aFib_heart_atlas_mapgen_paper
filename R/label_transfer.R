library(Seurat)

options(future.globals.maxSize = 8000 * 1024^2)

dat_new <- readRDS("/project2/xinhe/xsun/heart_atlas/2.label_transfer/data/global_raw.rds")
dat_new <- FindVariableFeatures(dat_new, selection.method = "vst", nfeatures = 2000, verbose = T)
saveRDS(dat_new,file = "/project2/xinhe/xsun/heart_atlas/2.label_transfer/data/dat_new_add_variable2.rds")

dat_our <- readRDS("/project2/gca/aselewa/heart_atlas_project/seurat/Heart_RNA_Processed_Combined_NoAtrium.rds")
dat_our$celltypes <- Idents(dat_our)

dat.anchors <- FindTransferAnchors(reference = dat_new, query = dat_our,dims = 1:30)
saveRDS(dat.anchors,file = "/project2/xinhe/xsun/heart_atlas/2.label_transfer/data/dat_anchor_transfer2.rds")

predictions <- TransferData(anchorset = dat.anchors, refdata = dat_new$cell_type, dims = 1:30)
saveRDS(predictions,file = "/project2/xinhe/xsun/heart_atlas/2.label_transfer/data/predictions2.rds")

dat_our <- AddMetaData(dat_our, metadata = predictions)
dat_our$prediction.match <- dat_our$predicted.id == dat_our$seurat_clusters
table(dat_our$prediction.match)
saveRDS(dat_our,file = "/project2/xinhe/xsun/heart_atlas/2.label_transfer/data/dat_our_after_transfer2.rds")



# Performing PCA on the provided reference using 1714 features as input.
# Projecting cell embeddings
# Finding neighborhoods
# Finding anchors
# Found 47507 anchors
# Filtering anchors
# Retained 6237 anchors