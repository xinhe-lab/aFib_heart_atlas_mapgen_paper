dat_new <- readRDS("/project2/xinhe/xsun/heart_atlas/2.label_transfer/data/dat_new_add_variable.rds")
dat_our <- readRDS("/project2/gca/aselewa/heart_atlas_project/seurat/Heart_RNA_Processed_Combined_NoAtrium.rds")
dat_our$celltypes <- Idents(dat_our)
prediction<- readRDS("/project2/xinhe/xsun/heart_atlas/2.label_transfer/data/prediction.rds")

celltypes_origin <- dat_our$celltypes

summary <- data.frame(cbind(names(celltypes_origin),as.character(celltypes_origin)))
colnames(summary) <- c("cell","celltype_our")

celltypes_predict <- data.frame(rownames(prediction),prediction$predicted.id)
summary <- merge(summary,celltypes_predict,by.x = "cell", by.y = "rownames.prediction.")
colnames(summary)[3] <- c("celltype_pred")

summary$celltype_match <- summary$celltype_pred

summary$celltype_match[summary$celltype_match == "Smooth_muscle_cells"] <- "Smooth Muscle"
summary$celltype_match[summary$celltype_match == "Pericytes"] <- "Pericyte"
summary$celltype_match[summary$celltype_match == "Ventricular_Cardiomyocyte"] <- "Cardiomyocyte"

summary$predict_match <- summary$celltype_our==summary$celltype_match
table(summary$predict_match)
# FALSE  TRUE 
# 1324 48035

index <- which(summary$celltype_match == "NotAssigned")
length(index)
# 349 

celltypes_ours <- levels(as.factor(summary$celltype_our))


#####heatmap

