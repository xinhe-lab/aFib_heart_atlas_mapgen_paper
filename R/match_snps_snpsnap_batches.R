#!/usr/bin/env Rscript

## match control SNPs using SNPsnap database
library(optparse)
library(foreach)
library(doParallel)
# library(data.table)

option_list = list(
  make_option("--input_file", action="store", default=NA, type='character',
              help="Prefix of PLINK files for the genotype [required]"),
  make_option("--dir_out", action="store", default=NA, type='character',
              help="Directory of output [required]"),
  make_option("--population", action="store", default="EUR", type='character',
              help="Population [optional]"),
  make_option("--num_matching_snps", action="store", default=100, type='integer',
              help="Number of matching SNPs [optional]"),
  make_option("--thresh_MAF", action="store", default=0.05, type='double',
              help="Minor Allele Frequency [optional]"),
  make_option("--thresh_ld_friends", action="store", default=0.5, type='double',
              help="LD buddies [optional]"),
  make_option("--num_batches", action="store", default=1, type='integer',
              help="Number of batches [optional]"),
  make_option("--idx_batch", action="store", default=1, type='integer',
              help="index of batches [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

# print(opt)

filename_input <- opt$input_file
dir_out <- opt$dir_out
population <- opt$population
num_matching_snps <- opt$num_matching_snps
thresh_MAF <- opt$thresh_MAF
thresh_ld_friends <- opt$thresh_ld_friends
num_batches <- opt$num_batches
idx_batch <- opt$idx_batch

##### Functions #####
## match control SNPs using SNPsnap database
## In five uniformly spaced increments, increase the allowable deviation for each of the genetic properties, ending with the prespecified maximum allowable deviation.
## For each increment, identify matching SNPs, defined as SNPs with genetic properties within the allowable deviations.
get_matching_SNPs <- function(snpsnap_input, snpsnap_db, exclude_HLA = TRUE, thresh_MAF, thresh_dist2peak, thresh_ld_friends){
  maf_range <- c(snpsnap_input$snp_maf - thresh_MAF, snpsnap_input$snp_maf + thresh_MAF)
  ld_friends_range <- snpsnap_input$friends_ld05 * c(1-thresh_ld_friends, 1+thresh_ld_friends)

  ## meet SNP Matching Criteria or not
  match_maf <- with(snpsnap_db, snp_maf >= maf_range[1] & snp_maf <= maf_range[2])
  match_ld <- with(snpsnap_db, friends_ld05 >= ld_friends_range[1] & friends_ld05 <= ld_friends_range[2])

  ## get matched indices
  idx_matched <-  which(match_maf & match_ld)
  snpID_matched <- snpsnap_db[idx_matched, snpID_type]

  ## exclude all the input SNPs from matching SNPs
  snpID_matched <- setdiff(snpID_matched, snpID_list_input_exist)

  ## exclude HLA SNPs (chr6 25Mb to 35Mb)
  if(exclude_HLA){
    snpID_inHLA <- snpsnap_db[snpsnap_db$inHLA == TRUE, snpID_type]
    snpID_matched <- setdiff(snpID_matched, snpID_inHLA)
  }

  return(snpID_matched)
}

match_snps_db <- function(snpsnap_input, snpID_list_input_exist, exclude_HLA = TRUE){

  k = 1
  while(k < 6){
    snpID_matched <- get_matching_SNPs(snpsnap_input, snpsnap_db, exclude_HLA, thresh_MAF_array[k], thresh_dist2peak_array[k], thresh_ld_friends_array[k])
    # cat("k=", k,":", length(snpID_matched), "SNPs matched \r")
    if(length(snpID_matched) >= num_matching_snps){
      break
    }else{
      k <- k + 1
    }
  }
  # cat("matching SNPs with", k, "increments \n")

  if(length(snpID_matched) >= num_matching_snps){
    ## matched enough control SNPs
    # cat("matched enough control SNPs \n")
    snpID_matched_selected <- sample(snpID_matched, num_matching_snps, replace = FALSE)
  }else if(length(snpID_matched) > 0 && length(snpID_matched) < num_matching_snps){
    ## not enough control SNPs
    # cat("not enough control SNPs matched\n")
    snpID_matched_selected <- sample(snpID_matched, num_matching_snps, replace = TRUE)
  }else{
    # cat("no control SNPs matched \n")
    snpID_matched_selected <- rep(NA, num_matching_snps)
  }

  return(c(snpsnap_input$snpID, length(snpID_matched), snpID_matched_selected))

}

split_SNPs_batches <- function(filename_input, dir_out, num_batches = 10){

  chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))

  prefix_name <- tools::file_path_sans_ext(basename(filename_input))
  dir_output <- paste0(dir_out, "/matched_SNPs/", prefix_name, "_", num_batches, "batches")
  dir.create(dir_output, showWarnings = F, recursive = T)

  snpID_list_input <- read.table(filename_input, header = F, stringsAsFactors = F)[,1]
  cat(length(snpID_list_input), "input SNPs. \n")

  snpID_list_input <- gsub("chr", "", snpID_list_input)

  if(any(grep("^rs", snpID_list_input))){
    snpID_type <- "rsID"
  }else if(any(grep("^[0-9]+:[0-9]+", snpID_list_input))){
    snpID_type <- "snpID"
  }else{
    print(head(snpID_list_input))
    stop("Please use correct snpID or rsID format! \n")
  }

  if(length(snpID_list_input) < num_batches){
    stop("Number of SNPs < number of batches!")
  }

  if(num_batches > 1){
    snpID_batches.l <- chunk(snpID_list_input, num_batches)
  }else{
    snpID_batches.l <- list(snpID_list_input)
  }

  saveRDS(snpID_batches.l, paste0(dir_output, "/", prefix_name, "_", num_batches, "batches.rds"))

  cat("\nSplit SNPs into", num_batches, "batches:", dir_output, "\n")

}

##### parameters #####
# directory for SNPsnap database
dir_SNPsnap_database <- "/project2/xinhe/kevinluo/SNPsnap/database/"

# use the environment variable SLURM_NTASKS_PER_NODE to set the number of cores
num_cores <- Sys.getenv("SLURM_NTASKS_PER_NODE")

if(is.na(num_cores) || num_cores == ""){
  num_cores <- 4
}

registerDoParallel(cores=num_cores)

cat("Using parallel", getDoParWorkers(), "cores \n")

prefix_name <- tools::file_path_sans_ext(basename(filename_input))
dir_output <- paste0(dir_out, "/matched_SNPs/", prefix_name, "_", num_batches, "batches")

cat("Input file:", filename_input, "\n")
cat("Population:", population, "\n")
cat("Find:", num_matching_snps, "matching SNPs", "\n")
cat("Threshold:",
    "\n thresh_MAF: ", thresh_MAF,
    "\n thresh_ld_friends: ", thresh_ld_friends, "\n")
cat("Batch #", idx_batch, "out of", num_batches, "\n")
cat("Output directory:", dir_output, "\n")

thresh_MAF_array <- seq(thresh_MAF/5, thresh_MAF, length.out = 5)
thresh_ld_friends_array <- seq(thresh_ld_friends/5, thresh_ld_friends, length.out = 5)

##### read input SNPs #####
filename_snp_batches <- paste0(dir_output, "/", prefix_name, "_", num_batches, "batches.rds")

if(!file.exists(filename_snp_batches)){
  cat(filename_snp_batches, "does not exist! \n")
  split_SNPs_batches(filename_input, dir_out, num_batches)
}

snpID_batches.l <- readRDS(filename_snp_batches)
snpID_list_input <- as.character(unlist(snpID_batches.l))
snpID_list_input_batch <- as.character(snpID_batches.l[[idx_batch]])

snpID_list_input <- gsub("chr", "", snpID_list_input)
snpID_list_input_batch <- gsub("chr", "", snpID_list_input_batch)

# head(snpID_list_input)
cat(length(snpID_list_input), "total input SNPs. \n")

if(any(grep("^rs", snpID_list_input))){
  snpID_type <- "rsID"
}else if(any(grep("^[0-9]+:[0-9]+", snpID_list_input))){
  snpID_type <- "snpID"
}else{
  stop("Please use correct snpID or rsID format! \n")
}

##### Begin matching SNPs #####
cat("Match control SNPs using SNPsnap database \n")

cat("Load SNPsnap database ... \n")

# snpsnap_db <- as.data.frame(fread(paste0(dir_SNPsnap_database, "/", population, "/ld0.5_collection.tab.gz")))
if(!file.exists(paste0(dir_SNPsnap_database, "/", population, "/ld0.5_collection.rds"))){
  stop("SNPsnap database does not exist!")
}

snpsnap_db <- readRDS(paste0(dir_SNPsnap_database, "/", population, "/ld0.5_collection.rds"))
rownames(snpsnap_db) <- snpsnap_db[,snpID_type]

if(! "inHLA" %in% colnames(snpsnap_db)){
  snpsnap_db_chr <- sapply(strsplit(snpsnap_db$snpID, ":"), "[[", 1)
  snpsnap_db_pos <- sapply(strsplit(snpsnap_db$snpID, ":"), "[[", 2)
  ## HLA region: chr6 25Mb - 35 Mb
  snpsnap_db$inHLA <- (as.character(snpsnap_db_chr) == "6" & (as.integer(snpsnap_db_pos) > 25e6 & as.integer(snpsnap_db_pos) < 35e6))
}

snpID_list_input_exist <- intersect(snpID_list_input, snpsnap_db[, snpID_type])

snpID_list_input_nonexist <- setdiff(snpID_list_input, snpsnap_db[, snpID_type])

cat(length(snpID_list_input_exist), "SNPs included in database. \n")

write.table(snpID_list_input_nonexist, paste0(dir_output, "/", "unincluded_snps.txt"), col.names = F, row.names = F, quote = F)

snpsnap_anno_input_exist <- snpsnap_db[snpID_list_input_exist, ]
write.table(snpsnap_anno_input_exist, paste0(dir_output, "/", "snpsnap_input_exist.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

cat("Matching SNPs for batch", idx_batch, "... \n")

snpID_list_input_batch_exist <- intersect(snpID_list_input_batch, snpsnap_db[, snpID_type])
cat(length(snpID_list_input_batch_exist), "input SNPs in this batch included in database \n")

snpsnap_anno_input_exist <- snpsnap_db[snpID_list_input_batch_exist, ]

start.time <- proc.time()
matched_snps.df <- foreach(i = 1:nrow(snpsnap_anno_input_exist), .combine=rbind) %dopar% {
  match_snps_db(snpsnap_anno_input_exist[i,], snpID_list_input_exist, exclude_HLA = TRUE)
}

colnames(matched_snps.df) <- c("input_SNP", "num_matches", paste("Set",  1:num_matching_snps, sep = "_"))
end.time <- proc.time()

print(end.time - start.time)

write.table(matched_snps.df, paste0(dir_output, "/", "SNPmatch_", num_matching_snps, "_matchedSNPs", "_batch", idx_batch, ".txt"),
            col.names = F, row.names = F, quote = F, sep = "\t")

cat("\nResults saved at:", dir_output, "\n")




