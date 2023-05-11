library(finemappeR)
library(tidyverse)
library(gwasvcf)
library(VariantAnnotation)

args <- commandArgs(trailingOnly = T)
EUR_LD_1KG <- '/project2/xinhe/1kg/bigsnpr/EUR_variable_1kg.rds'

gwas.vcf <- readVcf(file = args[1])
sumstats <- vcf_to_tibble(gwas.vcf)

bigSNP <- bigsnpr::snp_attach(rdsfile = EUR_LD_1KG)

gwas <- RunCleaner(sumstats = sumstats, 
                   ColsToKeep = c('seqnames','start','ES','SE','REF','ALT','ID','LP'), 
                   bigSNP = bigSNP)

saveRDS(gwas, args[2])


