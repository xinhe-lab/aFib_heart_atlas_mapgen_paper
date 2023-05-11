library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(GenomicRanges)
library(liftOver)
library(dplyr)
library(Gviz)
library(tidyverse)

ideal_order <- c("Cardiomyocyte","Smooth Muscle","Pericyte","Endothelial","Fibroblast","Neuronal", "Lymphoid","Myeloid")
palette <- readRDS('/project2/gca/aselewa/heart_atlas_project/notebooks/palette.rds')

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

ggClean <- function(rotate_axis=FALSE){
  tm <- theme_bw() +
    theme(text = element_text(size=18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=1))
  if(rotate_axis){
    tm <- tm + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }

  return(tm)

}

LegendOff <- function(){
  theme(legend.position = "none")
}

qcVlnPlot <- function(df, x,y,fill){
  ggplot(df, aes_string(x=x, y=y, fill=fill)) + geom_violin() + geom_jitter(shape=16, size=0.5, position=position_jitter(0.4)) + ggClean(rotate_axis=T) + xlab("") + LegendOff()
}

pullGeneScoreMatrix <- function(archr_project){
  gene.score <- getMatrixFromProject(ArchRProj = archr_project, useMatrix = 'GeneScoreMatrix')
  gene.score.mat <- as.matrix(assays(gene.score)[[1]])
  row.names(gene.score.mat) <- rowData(gene.score)$name
  gene.score.mat.norm <- imputeMatrix(gene.score.mat, imputeWeights = getImputeWeights(archr_project))
  return(gene.score.mat.norm)
}

custom_archr_plot <- function(archr_project, reduction = 'UMAP', group.by="regions", alpha=1, pt.size=1, palette=NULL, legend=TRUE, label=FALSE){

  cellColData <- as.data.frame(archr_project@cellColData)

  if(reduction == 'UMAP'){
    reduc.df <- as.data.frame(archr_project@embeddings$UMAP$df)
  }
  if(reduction == 'TSNE'){
    reduc.df <- as.data.frame(archr_project@embeddings$TSNE$df)
  }
  colnames(reduc.df) <- c("x1","x2")

  reduc.df[,"meta"] <- cellColData[,group.by, drop=T]
  p <- ggplot(reduc.df, aes(x=x1, y=x2, color=meta)) + ggrastr::rasterise(geom_point(size=pt.size, alpha=alpha), dpi = 200) + theme_set(theme_grey()) + ggClean() +
    xlab(paste0(reduction,'1')) + ylab(paste0(reduction,'2'))

  if(!is.null(palette)){
    p <- p + scale_color_manual(values = palette)
  }

  if(!legend){
    p <- p + LegendOff()
  }

  if(label){
    text.pos <- suppressMessages(reduc.df %>% group_by(meta) %>% summarise(pos_x = mean(x1), pos_y = mean(x2)))
    p <- p + ggrepel::geom_text_repel(data = text.pos, aes(x = pos_x, y = pos_y, label = meta), color="black")
  }

  return(p)
}


custom_dim_plot <- function(seurat, reduction = 'UMAP', group.by="cellTypes", alpha=1, pt.size=1, palette=NULL, legend=TRUE, label=FALSE){

  cellColData <- as.data.frame(seurat@meta.data)

  if(reduction == 'UMAP'){
    reduc.df <- as.data.frame(seurat@reductions$umap@cell.embeddings)
  }
  colnames(reduc.df) <- c("x1","x2")

  reduc.df[,"meta"] <- cellColData[,group.by, drop=T]
  p <- ggplot(reduc.df, aes(x=x1, y=x2, color=meta)) + ggrastr::rasterise(geom_point(size=pt.size, alpha=alpha), dpi = 200) + theme_set(theme_grey()) + ggClean() +
    xlab(paste0(reduction,'1')) + ylab(paste0(reduction,'2'))

  if(!is.null(palette)){
    p <- p + scale_color_manual(values = palette)
  }

  if(!legend){
    p <- p + LegendOff()
  }

  if(label){
    text.pos <- suppressMessages(reduc.df %>% group_by(meta) %>% summarise(pos_x = mean(x1), pos_y = mean(x2)))
    p <- p + ggrepel::geom_text_repel(data = text.pos, aes(x = pos_x, y = pos_y, label = meta), color="black")
  }

  return(p)
}


getQCPlots <- function(archr_project){
  cellColData <- as.data.frame(archr_project@cellColData)
  p1 <- qcVlnPlot(df = cellColData, x = "regions", y = "nFrags", fill = "regions")
  p2 <- qcVlnPlot(df = cellColData, x = "regions", y = "BlacklistRatio", fill = "regions")
  p3 <- qcVlnPlot(df = cellColData, x = "regions", y = "TSSEnrichment", fill="regions")
  p4 <- qcVlnPlot(df = cellColData, x = "regions", y = "NucleosomeRatio", fill="regions")
  p <- p1 + p2 + p3 + p4 + plot_layout(nrow=1)
}


RenameIdentity <- function(idents, from, to){
  new.idents <- plyr::mapvalues(idents, from = from, to = to)
  return(new.idents)
}

GRToString <- function(gr){
  paste0(as.character(seqnames(gr)),':',start(gr),'-',end(gr))
}

StringToGR <- function(st){
  chr <- sub(pattern = ':.*', replacement = "", st)
  st.end <- sub(pattern = '.*:', replacement = "", st)
  s <- as.numeric(sub(pattern = '-.*', replacement = "", st.end))
  e <- as.numeric(sub(pattern = '.*-', replacement = "", st.end))
  gr <- GenomicRanges::GRanges(seqnames = chr, IRanges::IRanges(start = s, end = e))
  return(gr)
}

snpIDtoGR <- function(id){
  chr <- sub(pattern = '_.*', replacement = "", x = id)
  tmp <- sub(pattern = 'chr[0-9]+_', replacement = "", x = id)
  pos <- sub(pattern = '_.+', replacement = "", x = tmp)

  id.gr <- StringToGR(paste0(chr, ':',pos,'-',pos))

  return(id.gr)
}

removeOverlaps <- function(X, to.remove){
  hits <- GenomicRanges::findOverlaps(query = X, subject = to.remove)
  if(length(hits) > 0){
    X <- X[-queryHits(hits),]
  }
  X
}

subsetByOverlapProp <- function(q, s, minPoverlap, maxgap=0){

  hits <- GenomicRanges::findOverlaps(query = q, subject = s, maxgap = maxgap)
  overlaps <- pintersect(q[queryHits(hits)], s[subjectHits(hits)])
  percentOverlap <- width(overlaps) / width(q[queryHits(hits)])
  hits <- hits[percentOverlap >= minPoverlap]

  return(hits)
}

SNPGenomeDistrib <- function(snp.gr, genomic.annots) {

  annot <- plyranges::join_overlap_inner(snp.gr, genomic.annots)

  annot.freq <- as.data.frame(annot@elementMetadata) %>%
    group_by(SNP) %>%
    dplyr::count(type) %>%
    mutate(n_bin = 1*(n>0)) %>%
    group_by(SNP) %>%
    mutate(f = n_bin / sum(n_bin))

  snpsIn <- length(unique(annot$SNP))
  nIntergenic <- length(snp.gr) - snpsIn
  snp.dist <- c("intergenic"=nIntergenic)
  snp.dist["exon"] <- sum(annot.freq$f[annot.freq$type=="exon"])
  snp.dist["UTR"] <- sum(annot.freq$f[annot.freq$type=="UTR"])
  snp.dist["intronic"] <- sum(annot.freq$f[annot.freq$type=="intron"])
  snp.dist["promoter"] <- sum(annot.freq$f[annot.freq$type=="promoter"])

  snp.dist.df <- data.frame(freq=100*snp.dist/length(snp.gr), category=names(snp.dist))
  snp.dist.df$category <- unfactor(snp.dist.df$category)

  snp.dist.df

}

join_overlap_list <- function(gr.list, X){
  res.list <- list()
  for(n in names(gr.list)){
    res.list[[n]] <- plyranges::join_overlap_inner(X, gr.list[[n]])
  }
  return(res.list)
}


get_elements_overlap_snps <- function(snp.gr, annotations) {
    for (f in annotations) {
        name <- paste0(basename(f), "_d")
        curr <- rtracklayer::import(f, format = "bed")
        seqlevelsStyle(curr) <- "NCBI"
        curr <- GenomicRanges::reduce(curr)
        overlap.df <- plyranges::join_overlap_inner(curr, snp.gr) %>%
            as_tibble() %>% mutate(enhancer = paste0('chr',seqnames, ':', start, '-', end)) %>%
            dplyr::select(snp, enhancer)
        colnames(overlap.df) <- c("snp", sub('.bed','',basename(annotations)))
    }
    return(overlap.df)
}

get_snp_to_gene_dist <- function(snp.gr, gene.gr, genes=NULL){

    if(!is.null(genes)){
        gene.gr <- gene.gr[gene.gr$gene_name %in% genes,]

    }

}

get_inner_gr_distance <- function(x.gr){

    x.gr.split <- GenomicRanges::split(x = x.gr, seqnames(x.gr))
    chr.dist <- list()

    for(i in 1:length(x.gr.split)){
        curr.gr <- x.gr.split[[i]]
        midpoints <- (start(curr.gr) + end(curr.gr))/2
        dist.distribution <- rep(0, length(midpoints)-1)
        for(j in 2:length(midpoints)){
            dist.distribution[j-1] <- midpoints[j] - midpoints[j-1]
        }

        chr.dist[[i]] <- dist.distribution
    }

    return(unlist(chr.dist, use.names = F))
}

hg38ToHg19 <- function(gr){

    path <- system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
    ch <- import.chain(path)
    unlist(liftOver(gr, ch))

}


computeKNN <- function(
    data = NULL,
    query = NULL,
    k = 50,
    includeSelf = FALSE,
    ...
){
    if(is.null(query)){
        query <- data
        searchSelf <- TRUE
    }else{
        searchSelf <- FALSE
    }
    if(searchSelf & !includeSelf){
        knnIdx <- nabor::knn(data = data, query = query, k = k + 1, ...)$nn.idx
        knnIdx <- knnIdx[,-1,drop=FALSE]
    }else{
        knnIdx <- nabor::knn(data = data, query = query, k = k, ...)$nn.idx
    }
    knnIdx
}


source_overlapCNN <- function(){
    Rcpp::sourceCpp(file = 'R/KNN_Utils.cpp')
}


prop_overlap_links <- function(X_e, X_p, Y_e, Y_p){

    enhancer.hits <- GenomicRanges::findOverlaps(query = X_e, subject = Y_e, maxgap = 1500)
    promoter.hits <- GenomicRanges::findOverlaps(query = X_p, subject = Y_p, maxgap = 1500)

    ehits.str <- paste0(queryHits(enhancer.hits),'-',subjectHits(enhancer.hits))
    phits.str <- paste0(queryHits(promoter.hits),'-',subjectHits(promoter.hits))

    ehitsIn <- queryHits(enhancer.hits)
    ehitsIn <- ehitsIn[ehits.str %in% phits.str]
    ehitsIn <- unique(ehitsIn)

    length(ehitsIn)/length(X_e)
}


# only requires the region of interest
knownGeneObject <- function(curr.locus.gr, genome){

    ga.track <- GenomeAxisTrack()
    if(genome == "hg19"){
        txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    }
    if(genome == "hg38"){
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    }

    grtrack <- GeneRegionTrack(range = txdb,
                               genome = "hg19",
                               chromosome = as.character(seqnames(curr.locus.gr)),
                               start = start(curr.locus.gr),
                               end = end(curr.locus.gr), name = "Genes")

    symbol(grtrack) <- mapIds(org.Hs.eg.db::org.Hs.eg.db,
                              keys=sub("\\.\\d+$", "", gene(grtrack)),
                              keytype="ENTREZID", column="SYMBOL")

    symbol(grtrack) <- ifelse(is.na(symbol(grtrack)), gene(grtrack), symbol(grtrack))
    return(grtrack)
}


load_genomic_annots <- function(genomic.annots.file){

  if( !file.exists(genomic.annots.file) ){
    stop("Genomic annotation file is not availble!")
  }
  cat("Load genomic annotations ...\n")

  genomic.annots <- readRDS(genomic.annots.file)

  genomic.annots$promoters$tss <- NA
  plus_strand <- which(strand(genomic.annots$promoters) == "+")
  minus_strand <- which(strand(genomic.annots$promoters) == "-")
  genomic.annots$promoters$tss[plus_strand] <- end(genomic.annots$promoters[plus_strand])
  genomic.annots$promoters$tss[minus_strand] <- start(genomic.annots$promoters[minus_strand])

  return(genomic.annots)
}

load_pcHiC <- function(pcHiC.file){

  if( !file.exists(pcHiC.file) ){
    stop("pcHiC file is not availble!")
  }
  cat("Load pcHiC data ...\n")

  pcHiC <- fread(pcHiC.file)

  pcHiC <- pcHiC %>% select(Promoter, Interacting_fragment)
  # separte genes connecting to the same fragment
  pcHiC <- pcHiC %>% separate_rows(Promoter) %>% dplyr::rename(gene_name = Promoter)

  pcHiC <- pcHiC %>%
    separate(Interacting_fragment, c("otherEnd_chr", "otherEnd_start", "otherEnd_end"), "\\.") %>%
    mutate(otherEnd_start = as.numeric(otherEnd_start), otherEnd_end = as.numeric(otherEnd_end))

  pcHiC.gr <- makeGRangesFromDataFrame(pcHiC, seqnames.field = "otherEnd_chr", start.field = "otherEnd_start", end.field = "otherEnd_end", keep.extra.columns = TRUE)
  seqlevelsStyle(pcHiC.gr) <- "UCSC"

  return(pcHiC.gr)
}


load_ABC <- function(ABC.file, ABC.thresh = 0, full.element = TRUE, flank = 0){

  if( !file.exists(ABC.file) ){
    stop("ABC file is not availble!")
  }
  cat("Load ABC data ...\n")

  ABC <- fread(ABC.file)

  if(full.element){
    ABC <- ABC %>%
      separate(name, c(NA, "element_region"), sep = "\\|", remove = FALSE) %>%
      separate(element_region, c(NA, "element_location"), sep = "\\:") %>%
      separate(element_location, c("element_start", "element_end"), sep = "\\-") %>%
      mutate(start = as.numeric(element_start), end = as.numeric(element_end))
  }

  ABC <- ABC %>%
    dplyr::rename(gene_name = TargetGene) %>%
    filter(ABC.Score >= ABC.thresh)

  if(flank > 0){
    ABC$start <- ABC$start - flank
    ABC$end <- ABC$end + flank
  }

  # ABC <- ABC %>% select(chr, start, end, name, class, gene_name, TargetGeneTSS, ABC.Score)

  ABC.gr <- makeGRangesFromDataFrame(ABC, keep.extra.columns = TRUE)
  seqlevelsStyle(ABC.gr) <- "UCSC"

  return(ABC.gr)
}

load_narrowpeaks <- function(peaks.file){
  peaks <- fread(peaks.file)
  colnames(peaks) <- c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
  peaks.gr <- makeGRangesFromDataFrame(peaks, keep.extra.columns = TRUE)
  seqlevelsStyle(peaks.gr) <- "UCSC"

  return(peaks.gr)
}

# Load and process coaccessibility data
process_coaccess <- function(celltype = c("CM", "all"), cor_thresh = 0.7, dist_thresh = 1e6){

  celltype <- match.arg(celltype)

  cat("Load and process coaccessibility data: celltype =", celltype, ", cor_thresh =", cor_thresh, ", dist_thresh =", dist_thresh, "\n")
  coacc_dir <- "/project2/xinhe/kevinluo/gene-level-finemapping/aFib_gene_finemapping/coaccessibility"
  if(toupper(celltype) == "CM"){
    coaccess <- readRDS(file.path(coacc_dir, "Coaccess_CMs_Enhancer_Promoters_corr_cut_-1_maxDist_1Mb_hg19.gr.rds"))
  }else if(toupper(celltype) == "ALL"){
    coaccess <- readRDS(file.path(coacc_dir, "Coacc_AllCellTypes_Enhancers_Promoter_Pairs_corrMin_-1_maxDist_1Mb_hg19.gr.rds"))
  }

  coaccess$dist <- abs(round((start(coaccess)+end(coaccess))/2) - round((coaccess$promoter_start + coaccess$promoter_end)/2))
  coaccess <- coaccess[coaccess$correlation > cor_thresh & coaccess$dist < dist_thresh,]

  return(coaccess)
}
