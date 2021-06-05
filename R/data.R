#' Normalized RNA Expression of 10 TCGA bladder cancer tumors.
#'
#' A dataset containing the normalized (by variance stabilizing transformation,
#' vst - see references) RNAseq expression of 10 randomly sampled tumors from
#' the bladder cancer cohort of The Cancer Genome Atlas.
#'
#' To reduce size of dataset, the genes that expressed in the lowest 20% (taking
#' ALL tumors into account, before random sampling) were removed from the
#' dataset.
#'
#' @format A data frame with 30773 rows and 10 columns
#' @source \url{https://portal.gdc.cancer.gov/projects/TCGA-BLCA}
#' @references Simon Anders, Wolfgang Huber: Differential expression analysis
#'   for sequence count data. Genome Biology 2010, 11:106.
#'   http://dx.doi.org/10.1186/gb-2010-11-10-r106
"blca"