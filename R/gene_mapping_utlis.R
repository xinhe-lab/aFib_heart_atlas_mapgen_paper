
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

#' Find all genes (promoters) within 1MB genes and assign weights based on distance
#'
#' @param snp.gr Genomic Ranges with the SNP locations.
#' @param promoters.gr Genomic Ranges with the gene promoter locations.
#' @param c.dist A scaling number used for computing weight based on SNP-gene distance. Weight = exp(-dist/c). Default = 1e5 (100kb).
#'
#' @return
gene_by_distance <- function(snp.gr, promoters.gr, c.dist = 1e5, dist.to = c("tss", "midpoint", "end")){

  res <- plyranges::join_overlap_inner(x = promoters.gr,
                                       y = snp.gr, maxgap = 1e6) # get all promoters within 1MB

  if(dist.to == "tss"){
    res$gene_pos <- res$tss
  }else if(dist.to == "midpoint"){
    res$gene_pos <- round((start(res) + end(res))/2)
  }else{
    res$gene_pos <- end(res)
  }

  res$distance <- abs(res$pos - res$gene_pos)
  res$weight <- exp(-res$distance / as.numeric(c.dist))

  return(res)
}


#' get gene credible sets from SNP level gene mapping table
#'
#' @param snp.gene.pip.mat A data frame of SNP level gene mapping table
#' @param by.locus If TRUE, get credible sets based on locus level gene PIP
#' @param gene.cs.percent.thresh percentage threshold for gene credible sets
#'
#' @return a data frame with gene credible sets
#' @export
get_gene_cs <- function(snp.gene.pip.mat, by.locus = TRUE, gene.cs.percent.thresh = 0.8){

  # add locus level gene PIP
  # For each locus - gene pair, sum over the fractional PIPs for SNPs in the locus and linked to the gene
  snp.locus.gene.pip.mat <- snp.gene.pip.mat %>%
    group_by(locus, gene_name) %>%
    mutate(locus_gene_pip = sum(pip * frac_pip)) %>% ungroup()

  # simplify to get locus, gene_name, locus_gene_pip, and gene_pip
  locus.gene.pip.df <- snp.locus.gene.pip.mat %>%
    dplyr::select(locus, gene_name, gene_pip, locus_gene_pip) %>%
    distinct(locus, gene_name, .keep_all=TRUE)

  # check if gene PIP is equal to the sum of gene-locus PIP
  for(gene in unique(locus.gene.pip.df$gene_name)){
    if(!all.equal(sum(locus.gene.pip.df$locus_gene_pip[locus.gene.pip.df$gene_name == gene]), locus.gene.pip.df$gene_pip[locus.gene.pip.df$gene_name == gene][1])){
      cat(gene, "'s gene_pip is not equal to the sum of locus_gene_pip! \n")
    }
  }

  if(by.locus){
    # for each locus, keep the genes with gene locus PIP cumsum > 0.8
    gene.cumsum.df <- locus.gene.pip.df %>%
      group_by(locus) %>%
      arrange(desc(locus_gene_pip)) %>%
      mutate(gene_pip_csum = cumsum(locus_gene_pip)) %>%
      dplyr::slice(1:which(gene_pip_csum >= gene.cs.percent.thresh)[1])

    # create gene cs table
    gene.cs.df <- gene.cumsum.df %>%
      group_by(locus) %>%
      summarise(gene_cs = paste0(gene_name, collapse=','),
                gene_cs_locus_pip = paste(paste0(gene_name, '(',round(locus_gene_pip,3),')'), collapse=','),
                top_gene = gene_name[1],
                top_locus_gene_pip = locus_gene_pip[1],
                top_gene_pip = gene_pip[1])
  }else{
    # for each locus, keep the genes with gene locus PIP cumsum > 0.8
    gene.cumsum.df <- locus.gene.pip.df %>%
      group_by(locus) %>%
      arrange(desc(gene_pip)) %>%
      mutate(gene_pip_csum = cumsum(gene_pip)) %>%
      dplyr::slice(1:which(gene_pip_csum >= gene.cs.percent.thresh)[1])

    # create gene cs table
    gene.cs.df <- gene.cumsum.df %>%
      group_by(locus) %>%
      summarise(gene_cs = paste0(gene_name, collapse=','),
                gene_cs_pip = paste(paste0(gene_name, '(',round(gene_pip,3),')'), collapse=','),
                top_gene = gene_name[1],
                top_gene_pip = gene_pip[1])
  }

  return(list(gene.cumsum.df = gene.cumsum.df, gene.cs.df = gene.cs.df, locus.gene.pip.df = locus.gene.pip.df))
}

