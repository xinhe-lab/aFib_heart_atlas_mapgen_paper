library(ArchR)
library(jsonlite)

ATAC_DIR <- '/project2/gca/Heart_Atlas/ATAC_seq/'
setwd('/project2/gca/aselewa/heart_atlas_project/ArchR/')
source('../R/analysis_utils.R')

# GLOBAL PARAMETERS
ATAC_SAMPLES <- c("MW200804AA","Deep_MW200804AC","Deep_MW200804AD","SP-HE-MW200928E1ATAC-413ATAC",
                 "SP-HE-HE200915ATAC-175ATAC","Deep_SP-HE200915ATAC-360ATAC","Deep_SP-HE200915ATAC-397ATAC","SP-HE-MW200928E1ATAC-398ATAC",
                 "SP-HE-MW200928E2ATAC-175ATAC","Deep_SP-MW200928E2ATAC-367ATAC","Deep_SP-MW200928E2ATAC-407ATAC","SP-HE-MW200928E1ATAC-408ATAC")

ATAC_INDIVIDUALS <- c(rep("02207",4),rep("02336",4),rep("03231",4))

ATAC_REGIONS <- rep(c("Septum","Right Ventricle","Left Ventricle","Apex"), 3)

# make arrow files

inputFiles <- sapply(ATAC_SAMPLES, function(x){paste0(ATAC_DIR,x,'/outs/fragments.tsv.gz')})
names(inputFiles) <- ATAC_SAMPLES

addArchRThreads(threads = 5)
addArchRGenome("hg38")
system('mkdir ArchR_ArrowFiles_raw_noAtrium')

setwd('ArchR_ArrowFiles_raw_noAtrium/')

for(i in 1:length(inputFiles)){
    ArrowFiles <- createArrowFiles(
        inputFiles = inputFiles[i],
        sampleNames = names(inputFiles)[i],
        outputNames = names(inputFiles)[i],
        minTSS = 1,
        minFrags = 100,
        verbose = T
    )
}

arrows <- list.files(path = '.', pattern = '*.arrow', full.names = F)
x <- addDoubletScores(input = arrows, k = 5)

setwd('../')


# create project

projDir <- 'ArchR_heart_noAtrium_raw'

projHeart <- ArchRProject(
    ArrowFiles = list.files(path = 'ArchR_ArrowFiles_raw_noAtrium', pattern = '*.arrow', full.names = T),
    outputDirectory = projDir,
    copyArrows = FALSE 
)
projHeart$individual <- plyr::mapvalues(x = projHeart$Sample, from = ATAC_SAMPLES, to = ATAC_INDIVIDUALS)
projHeart$regions <- plyr::mapvalues(x = projHeart$Sample, from = ATAC_SAMPLES, to = ATAC_REGIONS)

saveArchRProject(projHeart)

frag.sizes <- getFragmentSizes(projHeart)
saveRDS(frag.sizes, paste0(projDir, '/fragment_sizes.rds'))




