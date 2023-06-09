# Code and data resources for aFib heart atlas paper

This repository contains code and data resources to accompany our paper:
"Single-cell genomics improves the discovery of risk variants and genes of Atrial Fibrillation".

## Abstract

Genome-wide association studies (GWAS) have linked hundreds of loci to cardiac diseases. However, in most loci the causal variants and their target genes remain unknown. We developed a combined experimental and analytical approach that integrates single cell epigenomics with GWAS to prioritize risk variants and genes. We profiled accessible chromatin in single cells obtained from human hearts and leveraged the data to study genetics of Atrial Fibrillation (AF), the most common cardiac arrhythmia. Enrichment analysis of AF risk variants using cell-type-resolved open chromatin regions (OCRs) implicated cardiomyocytes as the main mediator of AF risk. We then performed statistical fine-mapping, leveraging the information in OCRs, and identified putative causal variants in 122 AF-associated loci. Taking advantage of the fine-mapping results, our novel statistical procedure for gene discovery prioritized 46 high-confidence risk genes, highlighting transcription factors and signal transduction pathways important for heart development. In summary, our analysis provides a comprehensive map of AF risk variants and genes, and a general framework to integrate single-cell genomics with genetic studies of complex traits.

## Reference

Alan Selewa\*, Kaixuan Luo\*, Michael Wasney, Linsin Smith, Xiaotong Sun, Chenwei Tang, Heather Eckart, Ivan Moskowitz, Anindita Basu, Xin He, Sebastian Pott. Single-cell genomics improves the discovery of risk variants and genes of Atrial Fibrillation.
medRxiv 2022.02.02.22270312; doi: https://doi.org/10.1101/2022.02.02.22270312

## Code and the Mapgen package

We implemented the computational workflow in the R package [mapgen][mapgen-link]. 
The version used in the paper: [doi:10.5281/zenodo.8067477][mapgen-zenodo-doi].

Additional analysis scripts and notebooks: 

 * [Analysis notebooks][Analysis-git-repo]
 * [R scripts][R-git-repo]

## Data

  * GEO repository for snRNA-seq and scATAC-seq data from our study: accession [GSE224997][GEO-link].

## License

All source code and software in this repository are made available under the terms of the MIT license.


[mapgen-link]: https://github.com/xinhe-lab/mapgen
[mapgen-zenodo-doi]: https://doi.org/10.5281/zenodo.8067477
[Analysis-git-repo]: https://github.com/xinhe-lab/aFib_heart_atlas_mapgen_paper/tree/main/analysis
[R-git-repo]: https://github.com/xinhe-lab/aFib_heart_atlas_mapgen_paper/tree/main/R
[GEO-link]: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224997
