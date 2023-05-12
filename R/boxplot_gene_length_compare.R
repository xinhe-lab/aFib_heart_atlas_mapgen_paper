library(biomaRt)
library(ggplot2)

setwd("/project2/xinhe/xsun/heart_atlas/1.OCR_compare")

human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

expr <- readRDS("~/seb/aselewa/heart_atlas_project/seurat/diff_expr_markers.df.rds")

celltypes <- c("Cardiomyocyte","Endothelial","Fibroblast")

gene_num <- 100
df_plot <- c()
for (i in 1:length(celltypes)){
  
  expr_cell <- expr[expr$cluster == celltypes[i],]
  expr_cell_sort <- expr_cell[order(expr_cell$avg_logFC,decreasing = T),]
  gene_list <- expr_cell_sort$gene[1:gene_num]
  
  gene_coords=getBM(attributes=c("hgnc_symbol","ensembl_gene_id", "start_position","end_position"),filters="hgnc_symbol", values=gene_list, mart=human)
  gene_length <- as.numeric(as.character(gene_coords$end_position)) -as.numeric(as.character(gene_coords$start_position)) +1
  
  df_tmp <- cbind(gene_coords,gene_length,rep(celltypes[i],nrow(gene_coords)))
  df_plot <- rbind(df_plot,df_tmp)
}
colnames(df_plot)[6] <- "celltypes"
df_plot$Study <- "Ours"

#save(df_plot,file = "df_genelength_100genes.rdata")





#human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

expr <- read.table("/project2/xinhe/xsun/heart_atlas/1.OCR_compare/expr_ren/Normalized.gene.expression.by.cluster.snRNA.txt",header = T)
celltypes <- c("Ventricular.Cardiomyocyte","Endothelial.cells","Fibroblast")
clusters<- c("vCardiomyocyte","Endothelial","Fibroblast")
expr_select <- expr[,colnames(expr) %in% celltypes]

gene_num <-100 
df_plot2 <- c()
for(i in 1:length(celltypes)) {
  
  expr_cell <- expr[,i]
  expr_cell <- data.frame(cbind(expr_cell,row.names(expr)))
  colnames(expr_cell) <- c("expr","gene") 
  
  expr_cell_sort <- expr_cell[order(as.numeric(as.character(expr_cell$expr)),decreasing = T),]
  gene_list <- expr_cell_sort$gene[1:gene_num]
  
  
  gene_coords=getBM(attributes=c("hgnc_symbol","ensembl_gene_id", "start_position","end_position"),filters="hgnc_symbol", values=gene_list, mart=human)
  gene_length <- as.numeric(as.character(gene_coords$end_position)) -as.numeric(as.character(gene_coords$start_position)) +1
  
  df_tmp <- cbind(gene_coords,gene_length,rep(clusters[i],nrow(gene_coords)))
  df_tmp <- merge(df_tmp,expr_cell_sort, by.x="hgnc_symbol",by.y="gene")
  df_plot2 <- rbind(df_plot2,df_tmp)
}

colnames(df_plot2)[6] <- "celltypes"
df_plot2$celltypes[df_plot2$celltypes == "vCardiomyocyte"] <-  "Cardiomyocyte"

df_plot2$celltypes <- factor(df_plot2$celltypes, levels=c("Cardiomyocyte","Endothelial","Fibroblast"))
df_plot2 <- df_plot2[,-7]
df_plot2$Study <- "Hocker's"

df_plot <- rbind(df_plot,df_plot2)

p <- ggplot(df_plot, aes(x=as.factor(celltypes), y=gene_length/1000, fill = Study)) + 
  geom_boxplot() +
  theme_bw(base_line_size =0.3) +
  labs(x = "Clusters", y= "Gene length (kb)") +
  theme(plot.title = element_text(size=35)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 35),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.title.y = element_text(size = 35),
        axis.text.y = element_text(size = 28, color = "black"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),) 

save(df_plot,file = "gene_length.rdata")
write.table(df_plot,file = "gene_length.txt", quote = F, row.names = F, col.names = T)
ggsave(filename = "plots/boxplot_genelenth_larger.pdf", device = "pdf", plot = p, dpi = 600,width = 10.5, height = 10, limitsize = FALSE )







