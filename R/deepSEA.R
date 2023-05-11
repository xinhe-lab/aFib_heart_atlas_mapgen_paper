

finemap.res <- read_csv('GWAS/gene_pips_summary_exp.csv')
finemap.res.topgene <- finemap.res %>% group_by(snp) %>% arrange(desc(gene_pip)) %>% slice(1) %>% filter(pip > 0.2) %>% arrange(desc(pip))

deepsea.res <- read_tsv('DeepSEA/Heart_Finemapped_SNPs/cb5429c9-6e7e-4b5b-ab27-99cdbed14583_high_conf_snps__1__FEATURE_evalue.tsv.gz')
snp.dis <- read_tsv('DeepSEA/Heart_Finemapped_SNPs/cb5429c9-6e7e-4b5b-ab27-99cdbed14583_high_conf_snps__1__VARIANT_dis.tsv')
snp.dis$snp <- deepsea.res$name

eval.thres <- 0.05
deepsea.evals <- as.data.frame(deepsea.res[,9:ncol(deepsea.res)])
rownames(deepsea.evals) <- deepsea.res$name
tissue.feature.names <- colnames(deepsea.evals)

relevant.cols <- c(contains(match = "ventricle", ignore.case = T, vars = tissue.feature.names),
                   contains(match = "atrium", ignore.case = T, vars = tissue.feature.names),
                   contains(match = "aorta", ignore.case = T, vars = tissue.feature.names),
                   contains(match = "heart", ignore.case = T, vars = tissue.feature.names),
                   contains(match = "muscle", ignore.case = T, vars = tissue.feature.names))

deepsea.evals.rel <- deepsea.evals[,relevant.cols]
deepsea.evals.melt <- reshape2::melt(t(deepsea.evals.rel))
colnames(deepsea.evals.melt) <- c("celltype_marker","snp","eval")
split.obj <- strsplit(as.character(deepsea.evals.melt$celltype_marker), split = '\\|')

deepsea.evals.melt$tissue_name <- sapply(split.obj, function(x){x[1]})
deepsea.evals.melt$marker_name <- sapply(split.obj, function(x){x[2]})

tissues.of.interest <- c("Left_Ventricle","Right_Ventricle","Right_Atrium","Aorta",
                         "Fetal_Heart","Muscle_Satellite_Cultured_Cells", "Skeletal_Muscle_Male","Skeletal_Muscle_Female",
                         "HSMM_Skeletal_Muscle_Myoblasts","HSMMtube_Skeletal_Muscle_Myotubes_Derived_from_HSMM")

deepsea.evals.melt <- deepsea.evals.melt[deepsea.evals.melt$tissue_name %in% tissues.of.interest,]

snp.tissue.summary <- deepsea.evals.melt %>% filter(eval <= eval.thres) %>% group_by(snp, tissue_name) %>% summarise(markers = paste0(marker_name, collapse='\n'))

snp.tissue.summary.wide <- reshape2::dcast(snp.tissue.summary, formula = snp ~ tissue_name)

final.res <- left_join(finemap.res.topgene, snp.tissue.summary.wide, on='snp') %>% 
    left_join(., snp.dis, on='snp') %>% 
    dplyr::select(-weight, -distance, -frac_pip, -fracExp, -chrom, -pos, -end) %>% 
    write_csv('GWAS/aFib_Finemapped_MinPip0.2_DeepSEA.csv')


