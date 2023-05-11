
setwd('/project2/gca/aselewa/heart_atlas_project/')
source('R/analysis_utils.R')

eqtls <- vroom::vroom('eQTL_enrich/summary_statistics/signifcant_summary_stats_v8/Heart_Left_Ventricle.v8.signif_variant_gene_pairs.txt.gz')

tmp <- sub(pattern = 'chr[0-9]+_', replacement = "", x = eqtls$variant_id)
pos <- sub(pattern = '_.+', replacement = "", x = tmp)
pos.int <- as.integer(pos)
bol <- !is.na(pos.int)
eqtls <- eqtls[bol,]
pos.int <- pos.int[bol]

chr <- sub(pattern = '_.*', replacement = "", x = eqtls$variant_id)
var.id.split <- strsplit(eqtls$variant_id, split = '_')
ref <- sapply(var.id.split, function(x){x[3]})
alt <- sapply(var.id.split, function(x){x[4]})
gene <- eqtls$gene_id
slope <- eqtls$slope
slope_se <- eqtls$slope_se
zscore <- eqtls$slope/eqtls$slope_se

eqtls.gr <- GRanges(seqnames = chr, ranges = IRanges(start = pos.int, end = pos.int), 
                    REF = ref, ALT = alt, 
                    gene = sub(pattern = '[.].*', '', gene),
                    slope = slope,
                    slope_se = slope_se,
                    zscore = zscore)

eqtls.gr.hg19 <- hg38ToHg19(eqtls.gr)

gene.annots <- readRDS('genomic_annotations/hg19_gtf_genomic_annots.gr.rds')$genes
gene.annots$gene_id_clean <- sub(pattern = '[.].*', '', gene.annots$gene_id)

eqtls.gr.hg19$Symbol <- gene.annots$gene_name[match(eqtls.gr.hg19$gene, gene.annots$gene_id_clean)]

seqlevelsStyle(eqtls.gr.hg19) <- "NCBI"

eqtls.gr.hg19$variant_ids <- paste0(as.character(seqnames(eqtls.gr.hg19)),'_',start(eqtls.gr.hg19),'_b37')

saveRDS(eqtls.gr.hg19, 'misc/V8_Signif_eQTLs_lifted_hg19.rds')



