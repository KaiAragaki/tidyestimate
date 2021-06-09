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
#' @param id either \code{entrezgene_id} or \code{hgnc_symbol}, whichever
#'   \code{df} contains.
#' @param tidy logical. If rownames contain gene identifier, set \code{FALSE}.
#'   If first column contains gene identifier, set \code{TRUE}
#' @param tell_missing logical. If \code{TRUE}, prints message of genes in common gene set that are not in supplied data frame.
#'
#' @return A data frame
#' @export
#'
#' @examples
#' filter_common_genes(ov, id = "hgnc_symbol", tidy = FALSE) |> 
#'   estimate_score(is_affymetrix = TRUE, tidy = FALSE)

filter_common_genes <- function (df, id = c("entrezgene_id", "hgnc_symbol"),
                                 tidy = FALSE, tell_missing = TRUE) {
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
        
        filtered
}
