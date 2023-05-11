
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
