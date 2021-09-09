#' Remove non-common genes from data frame
#'
#' As ESTIMATE score calculation is sensitive to the number of genes used, a set
#' of common genes used between six platforms has been established (see
#' \code{?tidyestimate::common_genes}). This function will filter for only those
#' genes.
#'
#' @param df a \code{data.frame} of RNA expression values, with columns corresponding
#'   to samples, and rows corresponding to genes. Either rownames or the first
#'   column can contain gene IDs (see \code{tidy})
#' @param id either \code{"entrezgene_id"} or \code{"hgnc_symbol"}, whichever
#'   \code{df} contains.
#' @param tidy logical. If rownames contain gene identifier, set \code{FALSE}.
#'   If first column contains gene identifier, set \code{TRUE}
#' @param tell_missing logical. If \code{TRUE}, prints message of genes in
#'   common gene set that are not in supplied data frame.
#' @param find_alias logical. If \code{TRUE} and \code{id = "hgnc_symbol"}, will
#'   attempt to find if genes missing from \code{common_genes} are going under
#'   an alias. See details for more information.
#'
#' @details
#'
#' The \code{find_aliases} argument will attempt to find aliases for HGNC
#' symbols in \code{tidyestimate::common_genes} but missing from the provided
#' dataset. This will only run if \code{find_aliases = TRUE} and \code{id =
#' "hgnc_symbol"}.
#'
#' This algorithm is very conservative: It will only make a match if the gene
#' from the common genes has only one alias that matches with only one gene from
#' the provided dataset, \emph{and} the gene from the provided dataset with
#' which it matches only matches with a single gene from the list of common
#' genes. (Note that a single gene may have many aliases). Once a match has been
#' made, the gene in the provided dataset is updated to the gene name in the
#' common gene list.
#'
#' While this method is fairly accurate, is is also a heuristic. Therefore, it is
#' disabled by default. Users should check which genes are becoming reassigned
#' to ensure accuracy.
#'
#' The method of generation of these aliases can be found at
#' \code{?tidyestimate::common_genes}
#'
#' @return A \code{tibble}, with gene identifiers as the first column
#' @export
#' @importFrom rlang .data
#' @examples
#' filter_common_genes(ov, id = "hgnc_symbol", tidy = FALSE, tell_missing = TRUE, find_alias = FALSE)

filter_common_genes <- function (df, id = c("entrezgene_id", "hgnc_symbol"),
                                 tidy = FALSE, tell_missing = TRUE, 
                                 find_alias = FALSE) {
  if (!tidy) {
    df <- dplyr::as_tibble(df, rownames = id)
  } 
  
  # Assumes first column contains gene IDs
  user_gene_col <- names(df)[1]

  if (sum(duplicated(df[[user_gene_col]])) > 0) warning("Input dataframe has duplicated IDs")
  
  common_genes <- tidyestimate::common_genes
  
  unique_common_genes <- unique(common_genes[,1:2])
  
  filtered <- dplyr::semi_join(df, common_genes, by = setNames("hgnc_symbol", user_gene_col)) 
  
  missing <- dplyr::anti_join(common_genes, df, by = setNames(user_gene_col, "hgnc_symbol"))[[user_gene_col]]
  
  if (find_alias & id == "hgnc_symbol") {
    filtered <- find_alias(df, common_genes, missing, filtered)
    missing <- dplyr::anti_join(common_genes, filtered, by = setNames(user_gene_col, "hgnc_symbol"))[["hgnc_symbol"]]
  }
  
  if (tell_missing) tell_missing(missing)
  
  message(glue::glue("
                     
                     Found {nrow(filtered)} of {nrow(unique_common_genes)} \\
                     genes ({round(nrow(filtered)/nrow(unique_common_genes) * 100, 3)}%) \\
                     in your dataset.
                     
                     ", sep = "\n\n"))
  
  filtered
}

tell_missing <- function(missing) {
  if (length(unique(missing)) > 10) {
    message(glue::glue("
                       
                       The following genes are in the list of common genes, \\
                       but not in your dataset:
                       
                       {paste(unique(missing)[1:10], collapse = ' ')}...and {length(unique(missing)) - 10} others
                       
                       ",
                       .sep = "\n\n"))
  } else{    
    message(glue::glue("
                       
                       The following genes are in the list of common genes, \\
                       but not in your dataset:
                       
                       {paste(unique(missing), collapse = ' ')}
                       
                       ",
                       .sep = "\n\n"))
  }
}

find_alias <- function(df, common_genes, missing, filtered) {

  df <- dplyr::mutate(df, row_id = rownames(df))

  matched_aliases <- common_genes |>
    dplyr::semi_join(as.data.frame(missing), 
                     by = c("hgnc_symbol" = "missing")) |>
    dplyr::inner_join(df, 
                      by = c("external_synonym" = "hgnc_symbol")) |> 
    dplyr::group_by(.data$hgnc_symbol) |>
    dplyr::mutate(hgnc_has_one_match = dplyr::if_else(dplyr::n() == 1, 
                                                            TRUE, 
                                                            FALSE)) |>
    dplyr::group_by(.data$row_id) |>
    dplyr::mutate(alias_has_one_match = dplyr::if_else(dplyr::n() == 1, 
                                                             TRUE, 
                                                             FALSE)) |>
    dplyr::filter(.data$hgnc_has_one_match & .data$alias_has_one_match) |>
    dplyr::rename(hgnc_symbol_new = .data$hgnc_symbol) |>
    dplyr::select(.data$row_id, .data$hgnc_symbol_new)
  
  df <- df |> 
    dplyr::right_join(matched_aliases, by = "row_id") |> 
    dplyr::mutate(alias = dplyr::if_else(!is.na(.data$hgnc_symbol_new), 
                                         .data$hgnc_symbol_new, 
                                         .data$hgnc_symbol)) |>
    dplyr::select(-.data$hgnc_symbol_new, -.data$hgnc_symbol, -.data$row_id) |>
    dplyr::relocate(.data$alias) |>
    dplyr::rename(hgnc_symbol = .data$alias)
  message(glue::glue("{nrow(df)} of {length(unique(missing))} missing genes found matches using aliases."))
  filtered <- dplyr::bind_rows(filtered, df)
}