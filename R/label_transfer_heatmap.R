library(ComplexHeatmap)

heat_map<-read.table("/project2/xinhe/xsun/heart_atlas/2.label_transfer/data/heatmap_matrix.txt",header = T)
#colnames(heat_map)[1] <- "Ventricular Cardiomyocyte"
colnames(heat_map)[1] <- "vCardiomyocyte"
colnames(heat_map)[8] <- "Smooth Muscle"

rownames(heat_map) <- c("Cardiomyocyte (8132)", "Endothelial (15336)","Fibroblast (11154)","Lymphoid (1002)","Myeloid (3356)","Neuronal (757)","Pericyte (7099)","Smooth Muscle (1376)")

pdf('/project2/xinhe/xsun/heart_atlas/2.label_transfer/heatmap.pdf', width=10, height=10)
#jpeg('/project2/xinhe/xsun/atlas/jcd_matrix_all_heatmap_uniq.jpg', width=1000, height=1000, res=200)
Heatmap(matrix = heat_map, cluster_rows = F, 
        cluster_columns = F, show_row_names = T, row_names_side = "left",
        name = "Proportion \nof matched",
        col = circlize::colorRamp2(c(0, 1), c("white","firebrick")), ###modify the figures
        row_title = NULL,
        column_title = NULL,
        row_gap = unit(1, "mm"),
        column_gap = unit(1, "mm"),
        na_col = "white",
        use_raster = T, 
        column_names_gp = grid::gpar(fontsize = 18),
        row_names_gp = grid::gpar(fontsize = 18))
dev.off()
