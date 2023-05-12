library(GenomicRanges)
library(rtracklayer)

setwd("/project2/xinhe/xsun/heart_atlas/5.OR")
load("/project2/xinhe/xsun/heart_atlas/5.OR/data/encode.beds.rdata")

### generate bed files, prepare for shuffling
# for (i in 1:length(encode.beds)){
#   
#   if (i <=6) {
#     gr <- encode.beds[[i]]
#   }else {
#     gr <- reduce(encode.beds[[i]])
#     gr[110637] <- shift(gr[110637],2)
#   }
#   
#   filename_bed <- paste0("data/random_prepare_samesize/",names(encode.beds)[i],"_random_prepare.bed")
#   export.bed(gr, con = filename_bed )
#   
# }


#### shuffling

file_prepare <- paste0("./data/random_prepare_samesize/",names(encode.beds),"_random_prepare.bed")
file_genome <- "/project2/xinhe/xsun/heart_atlas/3.ATAC_peak/random_ranges_data/my.genome"



for (i in 1:length(file_prepare)){
  
  for (j in 1:10) {
    
    file_exclude <- paste0("./data/gr_self/",names(encode.beds)[i],"_itself.bed" )
    file_out <- paste0("./data/data_shuffled_samesize_10times_exclude_itself/",names(encode.beds)[i],"_",j,"_shuffled.bed" )
    
    system(sprintf("/software/bedtools-2.27.1-el7-x86_64/bin/bedtools shuffle -i \"%s\" -g \"%s\" -excl \"%s\" -seed \"%s\" -noOverlapping |cat >\"%s\"",file_prepare[i], file_genome, file_exclude,j,file_out))
    
  }
  
}







