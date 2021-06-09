#' Remove non-common genes from data frame
#'
#' As ESTIMATE score calculation is sensitive to the number of genes used, a set
#' of common genes used between six platforms has been established (see
#' \code{?tidyestimate::common_genes} for more details). This function will
#' filter only for those genes.
#'
#' @param df a data frame of RNA expression, with columns corresponding to
#'   samples, and rows corresponding to genes. Either rownames or the first
#'   column can contain gene IDs (see \code{tidy})
#' @param id either \code{"entrezgene_id"} or \code{"hgnc_symbol"}, whichever
#'   \code{df} contains.
#' @param tidy logical. If rownames contain gene identifier, set \code{FALSE}.
#'   If first column contains gene identifier, set \code{TRUE}
#' @param tell_missing logical. If \code{TRUE}, prints message of genes in
#'   common gene set that are not in supplied data frame.
#' @param find_alias logical. If \code{TRUE}, will attempt to find if genes
#'   missing from \code{common_genes} are going under an alias. The function
#'   will attempt to find ESTIMATE gene set genes first (see
#'   \code{tidyestimate::gene_sets}), then move down the list to find any
#'   aliases for the remaining genes. See details for more information.
#'
#' @details
#'
#' The \code{find_aliases} argument will attempt to find aliases for HGNC
#' symbols in \code{tidyestimate::common_genes} but missing from the provided
#' dataset. This will only run if \code{find_aliases = TRUE} and \code{id =
#' "hgnc_symbol"}.
#'
#' Aliases will first attempt to be found for those in
#' \code{tidyestimate::gene_sets}, starting from the first missing gene in
#' \code{tidyestimate::common_genes} and moving downwards. The first encountered
#' match found in the given dataset with any alias of that given 'common gene'
#' will be marked as a hit, and the gene identifier of the dataframe provided
#' will be **changed** to the non-alias gene name in the common genes. The
#' matched alias will then be removed from all other aliases in order to
#' maximize the number of unique hits (as there are unique HGNC symbols that
#' share an alias). Once the process has been attempted for all missing genes of
#' the gene sets, attempts will be made for the remaining missing genes.
#'
#' As this can be error prone, it is disabled by default. Users should check
#' which genes are becoming reassigned to ensure accuracy.
#'
#' The method of generation of these aliases can be found at
#' \code{?tidyestimate::common_genes}
#'
#' @return A tidy tibble, with gene identifiers as the first column
#' @export
#'
#' @examples
#' filter_common_genes(ov, id = "hgnc_symbol", tidy, find_alias = FALSE) |>
#'   estimate_score(is_affymetrix = TRUE)

filter_common_genes <- function (df, id = c("entrezgene_id", "hgnc_symbol"),
                                 tidy = FALSE, tell_missing = TRUE, 
                                 find_alias = FALSE) {
        id <- match.arg(id)
        
        if (!tidy) {
                df <- dplyr::as_tibble(df, rownames = id)
        }
 
        common_genes <- tidyestimate::common_genes

        filtered <- dplyr::semi_join(df, common_genes) 

        message(glue::glue("Found {nrow(filtered)} of {nrow(common_genes)} genes ({round(nrow(filtered)/nrow(common_genes) * 100, 3)}%) in your dataset."))
        
        if (tell_missing) {
                missing <- dplyr::anti_join(common_genes, df)[[id]]
                message(glue::glue("The following genes are in the list of common genes, but not in your dataset:", 
                                   "{paste(missing, collapse = ' ')}", 
                                   .sep = "\n"))
        }
        
        if (find_alias) {
                # Take missing, then use aliases to search for hits. Start with
                # signature genes to make sure there's a higher chance that they
                # get matches. Then move on to other genes, going from top to
                # bottom. Be sure to eliminate matches that have already been
                # made as we go along.
                # Things to tell user:
                # Number of matches that have been made/salvaged using this method
                # Number of matches that were unable to be salvaged
                # Identity of matches that were unable to be salvaged
                # Then, write it all to the filtered tibble
                # Should only run if identifier is HGNC
                # Tell which genes are changed
        }
        
        filtered
}
