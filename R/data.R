#' Genes shared between six expression platforms
#'
#' The ESTIMATE model was trained on a set of genes shared between six
#' expression profiling platforms. Those genes are listed in this dataset.
#'
#' @description As the ESTIMATE model was trained on a specific set of genes,
#'   only those within this dataset should be included before running
#'   \code{estimate_scores}.
#'
#'   These are the genes common to 6 platforms:
#'
#'   - Affymetrix HG-U133Plus2.0
#'
#'   - Affymetrix HT-HG-U133A
#'
#'   - Affymetrix Human X3P
#'
#'   - Agilent 4x44K (G4112F)
#'
#'   - Agilent G4502A
#'
#'   - Illumina HiSeq RNA sequence
#'
#'   The Entrez IDs for the original 10412 genes were matched to HGNC symbols
#'   using \code{biomaRt}. Duplicates and blank entries were filtered. As some
#'   have now been discovered to be pseudogenes or have been deprecated, 22
#'   genes (at time of writing, June 2021) that were in the ESTIMATE package do
#'   not exist here.
#'
#'   As one gene can have multiple synonyms/aliases, and there is only one alias
#'   per line, the number of rows in the data frame (26339) does not reflect the
#'   number of unique genes in the dataset (10391).
#'
#' @format A data frame with 26339 rows and 3 variables: \describe{
#'   \item{entrezgene_id}{Entrez id of the gene} \item{hgnc_symbol}{Human Genome
#'   Organisation (HUGO) Gene Nomenclature Committee symbol}
#'   \item{external_synonym}{A synonym/alias a given gene may go by or
#'   previously went by} }
#' @source
#' \url{https://r-forge.r-project.org/scm/viewvc.php/pkg/estimate/data/common_genes.RData?root=estimate&view=log}
"common_genes"

#' Gene sets to infer tumor stromal and immune infiltration
#'
#' @description Two gene sets, each 141 genes in length, created to infer
#'   stromal and immune infiltration
#'   
#' @format A data frame with 141 row and 2 variables:
#'   \describe{
#'     \item{stromal_signature}{Geneset of HGNC symbols used to infer tumor stromal cell infiltration}
#'     \item{immune_signature}{Geneset of HGNC symbols used to infer tumor immune cell infiltration}
#'   }
#' @source 
#'   \url{https://r-forge.r-project.org/scm/viewvc.php/pkg/estimate/data/SI_geneset.RData?root=estimate&view=log}
"gene_sets"

#' Ovarian cancer tumor RNA expression
#'
#' @description A matrix containing RNA expression of 10 ovarian cancer tumors,
#'   measured using the Affymetrix U133Plus2.0 platform. These data have been
#'   rounded to the 4th decimal place to reduce file size.
#'
#' @format A matrix with 17256 rows and 10 columns, where each column represents
#'   a tumor, and each row represents a gene. Genes are represented by HGNC
#'   symbols in the rownames.
#' @source
#'   \url{https://r-forge.r-project.org/scm/viewvc.php/pkg/estimate/inst/extdata/sample_input.txt?root=estimate&view=log}
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
#' @source 
#'   \url{https://r-forge.r-project.org/scm/viewvc.php/pkg/estimate/data/PurityDataAffy.RData?root=estimate&view=log}
"purity_data_affy"
