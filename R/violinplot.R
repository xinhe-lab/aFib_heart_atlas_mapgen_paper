library(Seurat)
library(patchwork)
library(ggplot2)
library(ggrastr)

dat_our <- readRDS("/project2/gca/aselewa/heart_atlas_project/seurat/Heart_RNA_Processed_Combined_NoAtrium.rds")

# Relevel object@ident
dat_our_re <- dat_our
my_levels <- c('Cardiomyocyte', 'Endothelial','Fibroblast','Lymphoid','Myeloid','Neuronal','Pericyte','Smooth Muscle')
Idents(dat_our_re) <- factor(Idents(dat_our_re), levels= my_levels)

modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = FALSE,
                          plot.margin = margin(0, 0, 0, 0, "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )+
    ylab(feature) + 
    theme(legend.position = "none",
          plot.title= element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) + rasterise(geom_jitter(alpha=0.1, col="grey"),dpi = 100) ########rasterised
  p$layers <- p$layers %>% rev
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

#features<- c(  "TNNT2","VWF", "DCN","CD8A" ,"CD14", "PLP1","KCNJ8","MYH11")
features<- c(  "TNNT2","VWF", "DCN","CD3E" ,"FOLR2", "PLP1","KCNJ8","MYH11")
plot <-  StackedVlnPlot(obj = dat_our_re, features = features)

ggsave(filename = "violinplots_r50.pdf",plot = plot, dpi = 100, width = 10, height = 10)
