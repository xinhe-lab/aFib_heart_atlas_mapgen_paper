
## finemapping AF loci with TORUS prior using new ATAC peaks
library(mapgen)
library(dplyr)
suppressMessages(library(GenomicRanges))
suppressMessages(library(data.table))
suppressMessages(library(liftOver))
library(susieR)

setwd('/project2/xinhe/shared_data/aFib_gene_finemapping/pipeline/')

gwas.sumstats <- readRDS(paste0('GWAS/ebi-a-GCST006414_aFib.df.rds'))
head(gwas.sumstats)

# finemapping using SuSiE with L = 1 ----------------------------------------------------------------

# GWAS sample size
n = 1030836

# finemapping only significant GWAS Loci (P < 5e-8)
sig.loci <- gwas.sumstats %>%
  group_by(locus) %>%
  summarise(max_mlogP = max(pval)) %>%
  filter(max_mlogP > -log10(5e-8)) %>% pull(locus)

sumstats_finemap <- gwas.sumstats %>% dplyr::filter(locus %in% sig.loci) %>% dplyr::select(-torus_pip)

cat("Finemap",length(unique(sumstats_finemap$locus)), "loci.\n")

# Run SuSiE with L = 1
bigSNP <- bigsnpr::snp_attach(rdsfile = '/project2/xinhe/1kg/bigsnpr/EUR_variable_1kg.rds')

# susie_finemap_L1 is a list of SuSiE results, one for each chunk/LD block.
susie_finemap_L1 <- run_finemapping(sumstats_finemap, bigSNP, n = n, priortype = 'uniform', L = 1)

# add susie PIP information to GWAS summary stats
gwas_torus_finemapped <- merge_susie_sumstats(susie_results = susie_finemap_L1, sumstats = sumstats_finemap)
saveRDS(gwas_torus_finemapped, 'finemapping_latest_noAtrium_nonDA/AF_finemapping_result_unifprior_122loci.rds')

# Session Info -------------------
sessionInfo()

