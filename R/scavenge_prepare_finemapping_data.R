library(rtracklayer)
library(GenomicRanges)


#aFib_Finemapped_unif.tble <- readRDS("/project2/xinhe/shared_data/aFib_gene_finemapping/aFib_Finemapped_unif.tble.rds")
#aFib_Finemapped_unif.tble <- readRDS("/project2/xinhe/shared_data/aFib_gene_finemapping/aFib_Finemapped.tble.rds")
aFib_Finemapped_unif.tble <- readRDS("AF_finemapping_result_unifprior_122loci.rds")

finemap_hg19 <- data.frame(cbind(paste0("chr",aFib_Finemapped_unif.tble$chr),aFib_Finemapped_unif.tble$pos,aFib_Finemapped_unif.tble$pos,aFib_Finemapped_unif.tble$snp,aFib_Finemapped_unif.tble$susie_pip))
colnames(finemap_hg19) <- c("chr","start","end","snp","pip")

gr_hg19 <- GRanges(seqnames = finemap_hg19$chr, ranges = IRanges(names = finemap_hg19$snp, start = as.numeric(finemap_hg19$start), end = as.numeric(finemap_hg19$end)),score = finemap_hg19$pip)
chain <- import.chain("hg19ToHg38.over.chain")

result_hg38 <- as.data.frame(liftOver(gr_hg19,chain))
finemap_hg38 <- data.frame(cbind(result_hg38$seqnames,result_hg38$start,result_hg38$end,result_hg38$group_name,result_hg38$score))
finemap_hg38$X1 <- result_hg38$seqnames
write.table(finemap_hg38, file = "finemapping_aFib_hg38_new.bed", quote=F, col.names = F, row.names = F, sep="\t")
