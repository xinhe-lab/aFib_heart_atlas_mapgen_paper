library(GenomicFeatures)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
source('R/analysis_utils.R')

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(txdb) <- c(paste0("chr",1:22),"chrX")

entrez.to.symb <- function(IDs, db){
  AnnotationDbi::select(x = db, 
         keys = as.character(IDs), 
         columns = c("ENTREZID","SYMBOL"), 
         keytype = "ENTREZID")$SYMBOL
}

Tx.to.GeneID <- function(IDs, db){
  AnnotationDbi::select(x = db, 
         keys = as.character(IDs), 
         columns = "GENEID", 
         keytype = "TXID")$GENEID
}

add.tx.length <- function(IDs, db){
  res <- AnnotationDbi::select(x = db, 
         keys = as.character(IDs), 
         columns = c("TXSTART","TXEND"), 
         keytype = "TXID")
  res$TXEND - res$TXSTART
}

unlist_expand_ids <- function(gr){
  ids <- names(gr)
  expand.ids <- rep(ids, lengths(gr))
  gr <- unlist(gr)
  gr$tx_id <- expand.ids
  gr
}


annots <- list()

#genes
my.genes <- genes(txdb)
my.genes$gene_name <- entrez.to.symb(my.genes$gene_id, org.Hs.eg.db)

# get canonical transcript per gene
all.tx.ids <- AnnotationDbi::select(txdb, keys=my.genes$gene_id, columns=c("GENEID","TXID","TXSTART","TXEND"), keytype = "GENEID")
all.tx.ids$width <- all.tx.ids$TXEND - all.tx.ids$TXSTART
canonical_tx <- (all.tx.ids %>% group_by(GENEID) %>% summarise(canonical_txid = TXID[which.max(width)]))$canonical_txid

#exons
my.exons <- exonsBy(txdb, by = "tx")
my.exons <- unlist_expand_ids(my.exons)
my.exons <- my.exons[my.exons$tx_id %in% canonical_tx,]
my.exons$gene_id <- Tx.to.GeneID(my.exons$tx_id, db = txdb)
my.exons$gene_name <- entrez.to.symb(IDs = my.exons$gene_id, db = org.Hs.eg.db)

#promoters
my.promoters <- promoters(txdb, upstream = 2000, downstream = 0)
my.promoters <- my.promoters[my.promoters$tx_id %in% canonical_tx,]
my.promoters$gene_id <- AnnotationDbi::select(x = txdb, keys = as.character(my.promoters$tx_id), columns = "GENEID", keytype = "TXID")$GENEID
my.promoters$gene_name <- entrez.to.symb(my.promoters$gene_id, org.Hs.eg.db)

#introns (no need to annotate, we will grab nearest genes)
my.introns <- unlist_expand_ids(intronsByTranscript(txdb))
my.introns <- my.introns[my.introns$tx_id %in% canonical_tx,]
my.introns$gene_id <- Tx.to.GeneID(my.introns$tx_id, db = txdb)
my.introns$gene_name <- entrez.to.symb(IDs = my.introns$gene_id, db = org.Hs.eg.db)

# 5' UTRs
my.fiveUTRs <- unlist_expand_ids(fiveUTRsByTranscript(txdb))
my.fiveUTRs <- my.fiveUTRs[my.fiveUTRs$tx_id %in% canonical_tx,]
my.fiveUTRs$gene_id <- Tx.to.GeneID(my.fiveUTRs$tx_id, db = txdb)
my.fiveUTRs$gene_name <- entrez.to.symb(IDs = my.fiveUTRs$gene_id, db = org.Hs.eg.db)

# 3' UTRs
my.threeUTRs <- unlist_expand_ids(threeUTRsByTranscript(txdb))
my.threeUTRs <- my.threeUTRs[my.threeUTRs$tx_id %in% canonical_tx,]
my.threeUTRs$gene_id <- Tx.to.GeneID(my.threeUTRs$tx_id, db = txdb)
my.threeUTRs$gene_name <- entrez.to.symb(IDs = my.threeUTRs$gene_id, db = org.Hs.eg.db)

# add splice junctions
intron.chr <- c(as.character(seqnames(my.introns)), as.character(seqnames(my.introns)))
intron.pos <- c(start(my.introns), end(my.introns))
intron.gene.name <- c(my.introns$gene_name, my.introns$gene_name)
my.splice.junc <- GRanges(seqnames = intron.chr, ranges = IRanges(start = intron.pos-100, end = intron.pos+100))
my.splice.junc$gene_name <- intron.gene.name

annots <- list(exons = my.exons, 
               genes = my.genes,
               introns = my.introns,
               threeUTRs = my.threeUTRs,
               fiveUTRs = my.fiveUTRs,
               promoters = my.promoters,
               splice.junction = my.splice.junc)

# exons shouldnt contain UTRs
annots$exons <- removeOverlaps(X = annots$exons, to.remove = annots$threeUTRs)
annots$exons <- removeOverlaps(X = annots$exons, to.remove = annots$fiveUTRs)


saveRDS(annots, 'hg38_genomic_annotations.gr.rds')


