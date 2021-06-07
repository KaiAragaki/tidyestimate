#' Normalized RNA Expression of 10 TCGA bladder cancer tumors.
#'
#' A dataset containing the normalized (by variance stabilizing transformation,
#' `DESeq2::vst()`) RNAseq expression of 10 randomly sampled tumors from
#' the bladder cancer cohort of The Cancer Genome Atlas.
#'
#' To reduce size of dataset, the genes that expressed in the lowest 20% (taking
#' ALL tumors into account, before random sampling) were removed from the
#' dataset.
#'
#' @format A data frame with 30773 rows and 10 columns
#' @source \url{https://portal.gdc.cancer.gov/projects/TCGA-BLCA}
"blca"