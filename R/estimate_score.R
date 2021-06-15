#' Infer tumor purity using the ESTIMATE algorithm
#'
#' Infer tumor purity by using single-sample gene-set-enrichment-analysis with
#' stromal and immune cell signatures.
#'
#' @param df a \code{data.frame} of expression data, where columns are tumors
#'   and rows are genes. Gene names must be in the first column, and in the form
#'   of HGNC symbols.
#' @param is_affymetrix logical. Is the expression data from an Affymetrix
#'   array?
#'
#' @details
#'
#' ESTIMATE (and this tidy implementation) infers tumor infiltration using two
#' gene sets: a stromal signature, and an immune signature (see
#' \code{tidyestimate::gene_sets}).
#'
#' Enrichment scores for each sample are calculated using an implementation of
#' single sample Gene Set Enrichment Analysis (ssGSEA). Briefly, expression is
#' ranked on a per-sample basis, and the density and distribution of gene
#' signature 'hits' is determined. An enrichment of hits at the top of the
#' expression ranking confers a positive score, while an enrichment of hits at
#' the bottom of the expression ranking confers a negative score.
#'
#' An 'ESTIMATE' score is calculated by adding the stromal and immune scores
#' together.
#'
#' For Affymetrix arrays, an equation to convert an ESTIMATE score to a
#' prediction of tumor purity has been developed by Yoshihara et al. (see
#' references). It takes the approximate form of:
#'
#' \deqn{purity = cos(0.61 + 0.00015 * ESTIMATE)}
#'
#' Values have been rounded to two significant figures for display purposes.
#'
#'
#' @return A \code{data.frame} with sample names, as well as scores for stromal,
#'   immune, and ESTIMATE scores per tumor. If \code{is_affymetrix = TRUE},
#'   purity scores as well.
#' 
#' 
#' Purity scores can be interpreted absolutely: a purity of 0.9 means that tumor
#' is likely 90% pure by this metric. In the case where purity scores are not
#' available (such as in RNAseq), ESTIMATE scores can only be interpreted
#' relatively: a sample that has a lower ESTIMATE score than another in one
#' study can be regarded as more pure than another, but its absolute purity
#' cannot be inferred, nor can purity across other studies be inferred.
#'
#'
#' @references
#'
#' Barbie et al. (2009) <doi:10.1038/nature08460>
#'
#' Yoshihara et al. (2013) <doi:10.1038/ncomms3612>
#'
#'
#' @export
#' @importFrom rlang .data
#' @examples
#' filter_common_genes(ov, id = "hgnc_symbol", tidy = FALSE, tell_missing = TRUE, find_alias = TRUE) |> 
#'   estimate_score(is_affymetrix = TRUE)

estimate_score <- function(df, is_affymetrix) {

    rownames <- df[, 1]
    df <- df[, -1] 
    
    # Sample rank normalization
    df <- df |>  
        apply(2, as.numeric) |> 
        data.matrix()|> 
        apply(2, rank, ties.method = "average")
    
    rownames(df) <- unlist(rownames)
    
    df <- 10000*df/nrow(df)

    gene_set_names <- colnames(tidyestimate::gene_sets)
    scores <- matrix(NA_real_, nrow = ncol(df), ncol = 2)
    
    for (i in 1:2) {
        gene_set <- tidyestimate::gene_sets[, i, drop = FALSE]
        common_genes <- intersect(unlist(gene_set), rownames(df))
        
        message(glue::glue("
                           Number of {colnames(gene_set)} genes in data: \\
                           {length(common_genes)} (out of {nrow(gene_set)})
                           
                           ", .sep = "\n\n"))
        
        if (length(common_genes) == 0) { 
            next
        }
        
        ES.vector <- vector(length = ncol(df))
        
        for (sample in seq_along(1:ncol(df))) {
            
            # Arrange genes by INCREASING expression (remember, ranked previously)
            gene_list <- order(df[, sample], decreasing = TRUE)   
            ordered_genes <- df[gene_list, sample]
            
            hit_ind <- names(ordered_genes) %in% common_genes
            no_hit_ind <- 1 - hit_ind 
            
            ordered_genes <- ordered_genes^0.25
            
            hit_exp <- ordered_genes[hit_ind]
            
            no_hit_penalty <- 
                (no_hit_ind/sum(no_hit_ind)) |> 
                cumsum()
            
            hit_reward <- 
                ((hit_ind*ordered_genes)/sum(hit_exp)) |> 
                cumsum()

            ES.vector[sample] <- sum(hit_reward - no_hit_penalty)
        }
        scores[, i] <- ES.vector
    }

    scores <- scores |> 
        as.data.frame() |> 
        stats::setNames(c("stromal", "immune")) |> 
        dplyr::mutate(sample = colnames(df),
                      estimate = .data$stromal + .data$immune) |> 
        dplyr::relocate(sample)
   
    if (is_affymetrix) {
        # Calculate ESTIMATE-based tumor purity
        scores <- scores |> 
            dplyr::mutate(purity = cos(0.6049872018 + 0.0001467884 * .data$estimate),
                          purity = ifelse(.data$purity < 0, NA, .data$purity))
    }
    
    scores
}
