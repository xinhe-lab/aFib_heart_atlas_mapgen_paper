setwd('/project2/gca/aselewa/heart_atlas_project/')
source('R/analysis_utils.R')

union.ocrs <- readRDS('ArchR/ArchR_heart_latest_noAtrium/PeakCalls/UnionSet.gr.rds')
union.ocrs.hg19 <- hg38ToHg19(union.ocrs)

genomic.annots <- readRDS('genomic_annotations/hg19_gtf_genomic_annots.gr.rds')

# load HiC, prepare by running `prepare_HiC.R` script
genomic.annots$HiC <- readRDS('HiC/iPSC_CM_pcHiC_protein_Hg19.gr.rds') 

# load coaccessibility 
coaccess <- readRDS('ArchR/ArchR_heart_latest_noAtrium/CoAccessibility/Coaccess_CMs_corr_0_maxDist_1Mb_enhancers_hg19.gr.rds')
coaccess <- coaccess[coaccess$correlation>0.7,]
genomic.annots$coacc <- coaccess

# add active promoters
active.promoters <- subsetByOverlaps(genomic.annots$promoters, union.ocrs.hg19, minoverlap = 100)
genomic.annots$active.promoters <- active.promoters

# add nearby interactions
enhancers <- union.ocrs.hg19[union.ocrs.hg19$peakType!="Promoter",]

hits <- findOverlaps(enhancers, active.promoters, maxgap = 20000)
enhancer_nearby_promoter_20kb <- enhancers[queryHits(hits),]
enhancer_nearby_promoter_20kb$gene_name <- active.promoters$gene_name[subjectHits(hits)]

hits <- findOverlaps(enhancers, active.promoters, maxgap = 10000)
enhancer_nearby_promoter_10kb <- enhancers[queryHits(hits),]
enhancer_nearby_promoter_10kb$gene_name <- active.promoters$gene_name[subjectHits(hits)]

genomic.annots$enhancer_nearby_promoter_20kb <- enhancer_nearby_promoter_20kb
genomic.annots$enhancer_nearby_promoter_10kb <- enhancer_nearby_promoter_10kb

# save
saveRDS(genomic.annots, 'genomic_annotations/GWAS_gene_mapping_annots_hg19.gr.rds')


