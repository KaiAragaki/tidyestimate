#' Genes shared between six expression platforms
#' 
#' @description 
#' As the ESTIMATE model was trained on a specific set of genes, only those
#' within this dataset should be included before running \code{estimate_scores}.
#'
#' These are the genes common to 6 platforms:
#'
#' - Affymetrix HG-U133Plus2.0
#'
#' - Affymetrix HT-HG-U133A
#'
#' - Affymetrix Human X3P
#'
#' - Agilent 4x44K (G4112F)
#'
#' - Agilent G4502A
#'
#' - Illumina HiSeq RNA sequence
#'
#' The Entrez IDs for the original 10412 genes were matched to HGNC symbols
#' using \code{biomaRt}. Duplicates and blank entries were filtered. As some have now
#' been discovered to be pseudogenes or have been deprecated, 22 genes (at time
#' of writing, June 2021) that were in the ESTIMATE package do not exist here.
#'
#'
#' @format A data frame with 10391 rows and 2 variables: 
#' \describe{
#'   \item{entrezgene_id}{Entrez id of the gene} 
#'   \item{hgnc_symbol}{Human Genome Organisation (HUGO) Gene Nomenclature Committee symbol}
#'   }
#' @source \url{https://r-forge.r-project.org/scm/viewvc.php/pkg/estimate/data/common_genes.RData?root=estimate&view=log}
"common_genes"

#' Gene sets to measure tumor stromal and immune infiltration
#'
#' @description Two gene sets, each 141 genes in length, meant to measure
#'   stromal and immune infiltration
#'   
#' @format A data frame with 141 row and 2 variables:
#'   \describe{
#'     \item{stromal_signature}{Geneset of HGNC symbols used to measure stromal cell infiltration in a tumor}
#'     \item{immune_signature}{Geneset of HGNC symbols used to measure immune cell infiltration in a tumor}
#'   }
#' @source \url{https://r-forge.r-project.org/scm/viewvc.php/pkg/estimate/data/SI_geneset.RData?root=estimate&view=log}
"gene_sets"

#' Ovarian cancer tumor RNA expression
#'
#' @description A matrix containing RNA expression of 10 ovarian cancer tumors,
#'   measured using the Affymetrix U133Plus2.0 platform.
#'
#' @format A matrix with 17256 rows and 10 columns, where each column represents
#'   a tumor, and each row represents a gene. Genes are represented by HGNC
#'   symbols.
#' @source \url{https://r-forge.r-project.org/scm/viewvc.php/pkg/estimate/inst/extdata/sample_input.txt?root=estimate&view=log}
"ov"

#' Affymetrix data used to train ESTIMATE algorithm
#'
#' @description A data frame containing the ABSOLUTE-measured and
#'   ESTIMATE-predicted purity values of 995 tumors. Additionally, stromal and
#'   immune scores as calculated by ESTIMATE. All tumors were profiled on
#'   Affymetrix arrays, and were used to generate the Affymetrix algorithm.
#'
#' @format A data frame with 995 rows and 7 variables: 
#'   \describe{
#'     \item{purity_observed}{The purity of a tumor given by ABSOLUTE, ranging from 0 (least pure) to 1 (most pure)} 
#'     \item{stromal}{Stromal infiltration score, as measured by ESTIMATE}
#'     \item{immune}{Immune infiltration score, as measured by ESTIMATE}
#'     \item{estimate}{ESTIMATE score, calculated by the sum of immune and stromal scores}
#'     \item{purity_predicted}{Tumor purity inferred using the ESTIMATE algorithm}
#'     \item{ci_95_low}{Lower bound of a 95\% confidence interval of predicted purity scores}
#'     \item{ci_95_high}{Upper bound of a 95\% confidence interval of predicted purity scores}
#'   }
#'   
#' @source \url{https://r-forge.r-project.org/scm/viewvc.php/pkg/estimate/data/PurityDataAffy.RData?root=estimate&view=log}
"purity_data_affy"
