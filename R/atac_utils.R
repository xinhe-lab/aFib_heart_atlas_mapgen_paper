####################################################################
# Utils for generating Tile Matrix from 10X fragments file
####################################################################

library(tidyverse)
require(Rsamtools)
require(GenomicRanges)
require(S4Vectors)
require(BSgenome.Hsapiens.UCSC.hg38)

.getGenome <- function(genome){
  
  if(!(tolower(genome) %in% c('hg38'))){
    stop("Genome : ", genome, " is not currently supported.")  
  }
  if(tolower(genome)=="hg38"){
    if(!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)){
        stop("BSgenome for hg38 not installed! Please install by the following:\n\tBiocManager::install(\"BSgenome.Hsapiens.UCSC.hg38\")")
    }
    else{
      return(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
    }
  }
}

.chromLenGR <- function(genome, chrom){
  if(!(chrom %in% paste0('chr',1:22))){
    stop('Provide a valid chromosome: chr1-22')
  }
  genome <- .getGenome(genome)
  seqLen <- genome@seqinfo[chrom]@seqlengths
  chromGR <- GenomicRanges::GRanges(seqnames = chrom, IRanges::IRanges(start = 0, end = seqLen))
  
  return(chromGR)
}

.getFragments <- function(outs, params){
  tbx <- Rsamtools::TabixFile(paste0(outs,'fragments.tsv.gz'))
  tbx <- open.TabixFile(tbx)
  record <- Rsamtools::scanTabix(file = tbx, param = params)[[1]]

  if(length(record)>0){
    res <- record %>% 
      textConnection %>% 
      read.table(.) 
  }
  else{
    res <- NULL
  }
  return(res)
}

.getCellNames <- function(outs){
  singlecell <- read.table(paste0(outs,'singlecell.csv'), sep = ',', header = TRUE, stringsAsFactors = FALSE)
  barcodes <- singlecell[singlecell$is__cell_barcode==1,]$barcode
  return(barcodes)
}

.assignRegions <- function(chrom, binSize, nBins){
  GenomicRanges::GRanges(seqnames = rep(chrom, nBins),
                         IRanges::IRanges(start=seq(0, (nBins-1)*binSize, binSize), 
                                          end=seq(0, (nBins-1)*binSize, binSize)+binSize-1
                                          )
                         )
}

getTileMatrix10X <- function(outs, binSize=500, verbose=TRUE, binarize=TRUE, genome='hg38'){
  
  chroms <- sort(paste0("chr",1:22))
  
  singlecell <- read.table(paste0(outs,'singlecell.csv'), sep = ',', header = TRUE, stringsAsFactors = FALSE)
  singlecell <- singlecell[singlecell$is__cell_barcode==1,]
  cellNames <- singlecell$barcode
  
  matList <- list()
  regionList <- list()
  
  for(x in chroms){
    
    if(verbose){
      print(paste0('Getting Tile Matrix for ', x))
    }
    
    chromLen <- .chromLenGR(genome = 'hg38', chrom = x)
    nBins <- trunc(chromLen@ranges@width / binSize) + 1
    chromFrags <- .getFragments(outs, chromLen)
    bol <- chromFrags$V4 %in% cellNames
    if(sum(bol) > 0){
      
      regionList[[x]] <- .assignRegions(x, binSize, nBins)
      
      chromFrags <- chromFrags[chromFrags$V4 %in% cellNames,]
      cellIndex <- S4Vectors::match(chromFrags$V4, cellNames)
      
      mat <- Matrix::sparseMatrix(i = c(trunc(chromFrags$V2 / binSize), trunc(chromFrags$V3 / binSize)) + 1,
                                  j = c(cellIndex, cellIndex), 
                                  x = rep(1, 2*nrow(chromFrags)),
                                  dims = c(nBins, length(cellNames)))
      matList[[x]] <- mat
    }
  }

  tileMat <- Reduce(Matrix::rBind, matList)
  colnames(tileMat) <- cellNames
  regions <- GenomicRanges::bindROWS(x=NULL, regionList)
  rownames(tileMat) <- paste0(as.character(GenomicRanges::seqnames(regions)),':',GenomicRanges::start(regions),'-',GenomicRanges::end(regions))
  
  if(binarize){
    tileMat@x[tileMat@x > 0] <- 1
  }
  
  return(list(tileMat=tileMat, regions=regions, metadata=singlecell))
  
}




