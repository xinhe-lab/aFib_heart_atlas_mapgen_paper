setwd("/project2/xinhe/xsun/heart_atlas/5.OR/")

satac <- loadArchRProject('~/xsun_xin/heart_atlas/3.ATAC_peak/CopyOfArchR_Heart_Latest_Backup/')
peak.set <- satac@peakSet

markers <- readRDS('/project2/xinhe/xsun/heart_atlas/3.ATAC_peak/ArchR_Heart_Latest_Backup/PeakCalls/DA_MARKERS_FDRP_1_log2FC_1.rds')

#fnames <- c('ENCODE/H3K27ac/heart_left_ventricle.bed.gz', 'ENCODE/H3K27ac/heart_right_ventricle.bed.gz')
fnames <- c(list.files('data/ENCODE/DNase', pattern = "*.bed.gz$", full.names = T), 'data/ENCODE/H3K27ac/Heart_LVRV.bed.gz')
#'ENCODE/H3K27ac/heart_left_ventricle.bed.gz', 'ENCODE/H3K27ac/heart_right_ventricle.bed.gz')
encode.file.df <- data.frame(path = fnames, type = sub('_hg38','',sub('.bed.gz','',basename(fnames))), assay = basename(dirname(fnames)), stringsAsFactors = F)


encode.beds <- lapply(encode.file.df$path, function(x){ 
  encode <- vroom::vroom(x, col_names = F); 
  encode.gr <- GRanges(seqnames = encode$X1, ranges = IRanges(start = encode$X2, end = encode$X3));
  encode.gr
}
)
names(encode.beds) <- paste(encode.file.df$type)

encode.GR <- unlist(GRangesList(encode.beds))

markers$`Smooth Muscle` <- NULL
markers$Neuronal <- NULL

save(markers, file = "data/OCR_celltypespecific.rdata")
save(encode.GR, file = "data/encode.GR.rdata")
save(encode.file.df, file = "data/encode.file.df.rdata")
save(encode.beds, file = "data/encode.beds.rdata")
