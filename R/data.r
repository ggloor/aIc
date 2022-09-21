#' Selection-based differential sequence variant abundance dataset
#'
#' This data set gives the differential abundance of 1600 enzyme variants
#' grown under selective (NS) and selective (S) conditions
#'
#' @usage data(selex)
#' @format A data frame with 14 columns and 1600 rows
#' @source DOI:10.1073/pnas.1322352111
"selex"

#' Saccharomyces cerevisiae transcriptome
#'
#' A count table of a highly replicated RNA-seq experiment
#' with samples by column and genes by row. Two groups composed
#' of SNF2 knockout and WT, 48 samples in each.
#'
#' @usage data(transcriptome)
#' @format A data frame with 96 columns and 5892 rows
#' @source DOI: 10.1261/rna.053959.115 and PRJEB5348
"transcriptome"

#' single cell transcriptome data 
#'
#' A count table of a single cell transcriptome data subset from 
#' the count table from  doi:10.1038/s41592-019-0372-4. Two groups
#' memory T cells, and cytotoxic T cells, 1000 cells per group.
#' samples are by column and genes are by row. 
#'
#' @usage data(singleCell)
#' @format A data frame with 2000 columns and 1508 rows
#' @source https://www.nature.com/articles/s41592-019-0372-4
"singleCell"

#' meta-transcriptome data 
#'
#' A count table of a mixed population or metatranscriptome 
#' experiment. Two groups, H and BV are represented with 7 and 10
#' samples per group respectively. samples are by column and 
#' functions are by row.  
#'
#' @usage data(metaTscome)
#' @format A data frame with 17 columns and 3647 rows
#' @source doi:10.1007/978-3-030-71175-7_17 and doi:10.1007/978-1-4939-8728-3_13
"metaTscome"

#' 16S rRNA tag-sequencing data 
#'
#' A count table of a 16S rRNA amplicon data
#' Two groups, pupils and centenarians are represented with 198
#' and 161 samples per group respectively. samples are by column  
#' and OTU ids are by row.  
#' @usage data(meta16S)
#' @format A data frame with 359 columns and 860 rows
#' @source doi: 10.1128/mSphere.00327-17
"meta16S"