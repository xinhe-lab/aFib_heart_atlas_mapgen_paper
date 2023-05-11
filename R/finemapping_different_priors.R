
## finemapping AF loci with different prior settings
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(GenomicRanges))
library(mapgen)

data.dir <- '/project2/xinhe/shared_data/aFib_gene_finemapping/pipeline'
bigSNP <- bigsnpr::snp_attach(rdsfile = '/project2/xinhe/1kg/bigsnpr/EUR_variable_1kg.rds')

gwas.sumstats <- readRDS(paste0(data.dir, '/GWAS/ebi-a-GCST006414_aFib.df.rds'))
head(gwas.sumstats)

# GWAS sample size
n = 1030836

sig.loci <- gwas.sumstats %>%
  group_by(locus) %>%
  summarise(max_mlogP = max(pval)) %>%
  filter(max_mlogP > -log10(5e-8)) %>% pull(locus)

gwas.sumstats.sigloci <- gwas.sumstats[gwas.sumstats$locus %in% sig.loci, ]
cat(length(unique(gwas.sumstats.sigloci$locus)), "GWAS significant loci. \n")

# save gene body bed annotation file
genomic.annots <- readRDS(file.path(data.dir, "annotations/GWAS_gene_mapping_annots_hg19.gr.rds"))
gene.annots <- genomic.annots$genes
gene.locations <- as.data.frame(gene.annots)[, c('seqnames', 'start', 'end', 'gene_name', 'strand')]
dir.create(paste0(data.dir, "/annotations/gene_annotations"), recursive = TRUE)
gene.locations.bed <- gene.locations[, c('seqnames', 'start', 'end', 'gene_name')] %>% dplyr::rename(chr = seqnames, name = gene_name)
fwrite(gene.locations.bed, file.path(data.dir, "annotations/gene_annotations/gene.locations.hg19.bed"), sep = "\t", col.names = FALSE)

gene.locations$tss <- GenomicRanges::start(GenomicRanges::resize(gene.annots, width = 1))
tss20kb.bed <- gene.locations
tss20kb.bed$start <- gene.locations$tss - 2e4
tss20kb.bed$end <- gene.locations$tss + 2e4
tss20kb.bed <- tss20kb.bed[, c('seqnames', 'start', 'end', 'gene_name')] %>% dplyr::rename(chr = seqnames, name = gene_name)
write.table(tss20kb.bed, file.path(data.dir, "annotations/gene_annotations/TSS.20kb.hg19.bed"), sep = "\t",
            col.names=FALSE, row.names=FALSE, quote=FALSE)


# TORUS prior -------------------------------------------------------------------------------------------
annotation_bed_files <- list.files(paste0(data.dir, '/finemapping_latest_noAtrium_nonDA/annotations_for_finemapping_hg19'), pattern = '*.bed', full.names = T)
torus.files <- prepare_torus_input_files(gwas.sumstats,
                                         annotation_bed_files,
                                         torus_input_dir = paste0(data.dir, '/finemapping_latest_noAtrium_nonDA/torus_input'))

torus.result <- run_torus(torus.files$torus_annot_file,
                          torus.files$torus_zscore_file,
                          option = "est-prior",
                          torus_path = "torus") # set the path to 'torus' executable.
torus.enrich <- torus.result$enrich
torus.prior <- torus.result$snp_prior
saveRDS(torus.result, paste0(data.dir, '/finemapping_latest_noAtrium_nonDA/Torus_Enrichment_Results_Joint.rds'))

sumstats.sigloci <- prepare_susie_data_with_torus_result(sumstats = gwas.sumstats.sigloci,
                                                         torus_prior = torus.prior)
sumstats.sigloci <- sumstats.sigloci %>% select(-torus_pip)
saveRDS(sumstats.sigloci, paste0(data.dir, '/finemapping_latest_noAtrium_nonDA/sumstats.sigloci.rds'))

sumstats.sigloci <- readRDS(paste0(data.dir, '/finemapping_latest_noAtrium_nonDA/sumstats.sigloci.rds'))
cat("Finemap",length(unique(sumstats.sigloci$locus)), "loci.\n")

# susie_finemap_L1 is a list of SuSiE results, one for each chunk/LD block.
susie_finemap_L1 <- run_finemapping(sumstats.sigloci, bigSNP, priortype = 'torus', L = 1)

# add susie PIP information to GWAS summary stats
gwas_torus_finemapped <- merge_susie_sumstats(susie_results = susie_finemap_L1, sumstats = sumstats.sigloci)
saveRDS(gwas_torus_finemapped, paste0(data.dir, '/finemapping_latest_noAtrium_nonDA/AF_finemapping_result_torusprior_122loci.rds'))

# Uniform prior -------------------------------------------------------------------------------------------
sumstats.sigloci <- readRDS(paste0(data.dir, '/finemapping_latest_noAtrium_nonDA/sumstats.sigloci.rds'))
cat("Finemap",length(unique(sumstats.sigloci$locus)), "loci.\n")

susie_finemap_L1 <- run_finemapping(sumstats.sigloci, bigSNP, priortype = 'uniform', L = 1)
gwas_unif_finemapped <- merge_susie_sumstats(susie_results = susie_finemap_L1, sumstats = sumstats.sigloci)
saveRDS(gwas_unif_finemapped, paste0(data.dir, '/finemapping_latest_noAtrium_nonDA/AF_finemapping_result_unifprior_122loci.rds'))


# Shuffled prior -------------------------------------------------------------------------------------------
sumstats.sigloci <- readRDS(paste0(data.dir, '/finemapping_latest_noAtrium_nonDA/sumstats.sigloci.rds'))
cat("Random shuffle torus prior ...\n")
set.seed(10)
sumstats.sigloci$torus_prior <- sample(sumstats.sigloci$torus_prior)

cat("Finemap",length(unique(sumstats.sigloci$locus)), "loci.\n")
susie_finemap_L1 <- run_finemapping(sumstats.sigloci, bigSNP, n = n, priortype = 'torus', L = 1)
gwas_shuffled_finemapped <- merge_susie_sumstats(susie_results = susie_finemap_L1, sumstats = sumstats.sigloci)
saveRDS(gwas_shuffled_finemapped, paste0(data.dir, '/finemapping_latest_noAtrium_nonDA/AF_finemapping_result_shuffledprior_122loci.rds'))

# Gene body set prior -------------------------------------------------------------------------------------------
sumstats.sigloci <- readRDS(paste0(data.dir, '/finemapping_latest_noAtrium_nonDA/sumstats.sigloci.rds'))

sumstats.sigloci.gr <- makeGRangesFromDataFrame(sumstats.sigloci, start.field = 'pos', end.field = 'pos')
seqlevelsStyle(sumstats.sigloci.gr) <- 'UCSC'

genomic.annots <- readRDS(file.path(data.dir, "annotations/GWAS_gene_mapping_annots_hg19.gr.rds"))
gene.annots <- genomic.annots$genes

in_genebody <- which(countOverlaps(sumstats.sigloci.gr, gene.annots, ignore.strand=TRUE) > 0)
cat(length(in_genebody), "snps in gene body. \n")

sumstats.sigloci$torus_prior <- 1e-5
sumstats.sigloci$torus_prior[in_genebody] <- 1e-4
range(sumstats.sigloci$torus_prior)

cat("Finemap",length(unique(sumstats.sigloci$locus)), "loci.\n")
susie_finemap_L1 <- run_finemapping(sumstats.sigloci, bigSNP, n = n, priortype = 'torus', L = 1)

gwas_genebodyprior_finemapped <- merge_susie_sumstats(susie_results = susie_finemap_L1, sumstats = sumstats.sigloci)
saveRDS(gwas_genebodyprior_finemapped, paste0(data.dir, '/finemapping_latest_noAtrium_nonDA/AF_finemapping_result_genebody.setprior_122loci.rds'))

# TSS 20kb set prior -------------------------------------------------------------------------------------------
sumstats.sigloci <- readRDS(paste0(data.dir, '/finemapping_latest_noAtrium_nonDA/sumstats.sigloci.rds'))

sumstats.sigloci.gr <- makeGRangesFromDataFrame(sumstats.sigloci, start.field = 'pos', end.field = 'pos')
seqlevelsStyle(sumstats.sigloci.gr) <- 'UCSC'

genomic.annots <- readRDS(file.path(data.dir, "annotations/GWAS_gene_mapping_annots_hg19.gr.rds"))
gene.annots <- genomic.annots$genes
gene.locations <- as.data.frame(gene.annots)[, c('seqnames', 'start', 'end', 'gene_name', 'strand')]
gene.locations$tss <- GenomicRanges::start(GenomicRanges::resize(gene.annots, width = 1))
tss20kb <- gene.locations
tss20kb$start <- gene.locations$tss - 2e4
tss20kb$end <- gene.locations$tss + 2e4
tss20kb.gr <- makeGRangesFromDataFrame(tss20kb)

in_TSS20kb <- which(countOverlaps(sumstats.sigloci.gr, tss20kb.gr, ignore.strand = TRUE) > 0)
cat(length(in_TSS20kb), "snps in 20kb from TSS. \n")

sumstats.sigloci$torus_prior <- 1e-5
sumstats.sigloci$torus_prior[in_TSS20kb] <- 1e-4
range(sumstats.sigloci$torus_prior)

cat("Finemap",length(unique(sumstats.sigloci$locus)), "loci.\n")
susie_finemap_L1 <- run_finemapping(sumstats.sigloci, bigSNP, n = n, priortype = 'torus', L = 1)

gwas_TSS20kb_setprior_finemapped <- merge_susie_sumstats(susie_results = susie_finemap_L1, sumstats = sumstats.sigloci)
saveRDS(gwas_TSS20kb_setprior_finemapped, paste0(data.dir, '/finemapping_latest_noAtrium_nonDA/AF_finemapping_result_TSS20kb.setprior_122loci.rds'))

# Gene body prior estimated by TORUS -----------------------------------------------------------------------------------
annotation_bed_files <- file.path(data.dir, "annotations/gene_annotations/gene.locations.hg19.bed")
torus.files <- prepare_torus_input_files(gwas.sumstats,
                                         annotation_bed_files,
                                         torus_input_dir = paste0(data.dir, '/finemapping_latest_noAtrium_nonDA/torus_input/genebody_est'))
torus_annot <- fread(paste0(data.dir, '/finemapping_latest_noAtrium_nonDA/torus_input/genebody_est/torus_annotations.txt.gz'))

range(torus_annot[,2])

torus.result <- run_torus(torus.files$torus_annot_file,
                          torus.files$torus_zscore_file,
                          option = "est-prior",
                          torus_path = "torus") # set the path to 'torus' executable.
torus.enrich <- torus.result$enrich
torus.prior <- torus.result$snp_prior
dir.create(paste0(data.dir, '/finemapping_latest_noAtrium_nonDA/torus_result/'))
saveRDS(torus.result, paste0(data.dir, '/finemapping_latest_noAtrium_nonDA/torus_result/torus_genebody_enrichment.rds'))

table(torus.prior$torus_prior)
torus.prior[torus.prior$snp == "rs115315451", ]

sumstats.sigloci <- prepare_susie_data_with_torus_result(sumstats = gwas.sumstats.sigloci,
                                                         torus_prior = torus.prior)

sumstats.sigloci <- sumstats.sigloci %>% select(-torus_pip)
saveRDS(sumstats.sigloci, paste0(data.dir, '/finemapping_latest_noAtrium_nonDA/sumstats.sigloci.genebody.est.rds'))

sumstats.sigloci <- readRDS(paste0(data.dir, '/finemapping_latest_noAtrium_nonDA/sumstats.sigloci.genebody.est.rds'))
cat("Finemap",length(unique(sumstats.sigloci$locus)), "loci.\n")

# susie_finemap_L1 is a list of SuSiE results, one for each chunk/LD block.
susie_finemap_L1 <- run_finemapping(sumstats.sigloci, bigSNP, n = n, priortype = 'torus', L = 1)

# add susie PIP information to GWAS summary stats
gwas_torus_finemapped <- merge_susie_sumstats(susie_results = susie_finemap_L1, sumstats = sumstats.sigloci)
saveRDS(gwas_torus_finemapped, paste0(data.dir, '/finemapping_latest_noAtrium_nonDA/AF_finemapping_result_genebody.est.prior_122loci.rds'))

# TSS 20kb prior estimated by TORUS -----------------------------------------------------------------------------------
annotation_bed_files <- file.path(data.dir, "annotations/gene_annotations/TSS.20kb.hg19.bed")
torus.files <- prepare_torus_input_files(gwas.sumstats,
                                         annotation_bed_files,
                                         torus_input_dir = paste0(data.dir, '/finemapping_latest_noAtrium_nonDA/torus_input/tss20kb_est'))

torus.result <- run_torus(torus.files$torus_annot_file,
                          torus.files$torus_zscore_file,
                          option = "est-prior",
                          torus_path = "torus") # set the path to 'torus' executable.
torus.enrich <- torus.result$enrich
torus.prior <- torus.result$snp_prior
saveRDS(torus.result, paste0(data.dir, '/finemapping_latest_noAtrium_nonDA/torus_result/torus_tss20kb_enrichment.rds'))

sumstats.sigloci <- prepare_susie_data_with_torus_result(sumstats = gwas.sumstats.sigloci,
                                                         torus_prior = torus.prior)

sumstats.sigloci <- sumstats.sigloci %>% select(-torus_pip)
saveRDS(sumstats.sigloci, paste0(data.dir, '/finemapping_latest_noAtrium_nonDA/sumstats.sigloci.genebody.est.rds'))

sumstats.sigloci <- readRDS(paste0(data.dir, '/finemapping_latest_noAtrium_nonDA/sumstats.sigloci.genebody.est.rds'))
cat("Finemap",length(unique(sumstats.sigloci$locus)), "loci.\n")

# susie_finemap_L1 is a list of SuSiE results, one for each chunk/LD block.
susie_finemap_L1 <- run_finemapping(sumstats.sigloci, bigSNP, n = n, priortype = 'torus', L = 1)
# add susie PIP information to GWAS summary stats
gwas_torus_finemapped <- merge_susie_sumstats(susie_results = susie_finemap_L1, sumstats = sumstats.sigloci)
saveRDS(gwas_torus_finemapped, paste0(data.dir, '/finemapping_latest_noAtrium_nonDA/AF_finemapping_result_tss20kb.est.prior_122loci.rds'))

