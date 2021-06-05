#' Calculate Stromal, Immune, and ESTIMATE scores
#'
#' @param df a dataframe of normalized expression data, where columns are
#'   patients and genes are rows. Gene names must be HGNC, and may either be the
#'   first column or the rownames (see `tidy` argument below)
#' @param is_affymetrix logical. Is the expression data from an Affymetrix
#'   array?
#' @param tidy logical. If true, assumes the first column contains HGNC symbols.
#'   If false, assumes rownames contain HGNC symbols
#'
#' @return A data frame with scores for stromal, immune, and ESTIMATE scores per tumor.
#' 
#' @export
#' @importFrom rlang .data

estimate_score <- function(df, is_affymetrix, tidy) {

    if(tidy) {
        rownames <- df[, 1]
        df <- df[, -1] 
    } else {
        rownames <- row.names(df)
    }
    
    # Sample rank normalization
    df <- df |>  
        apply(2, as.numeric) |> 
        data.matrix()|> 
        apply(2, rank, ties.method = "average")
    
    rownames(df) <- rownames
    
    df <- 10000*df/nrow(df)

    gene_set_names <- colnames(gene_sets)
    scores <- matrix(NA_real_, nrow = ncol(df), ncol = 2)
    
    for (i in 1:2) {
        gene_set <- gene_sets[, i, drop = FALSE]
        common_genes <- intersect(unlist(gene_set), rownames(df))
        
        message(glue::glue("Gene set: {colnames(gene_set)}",
                           "# gene set genes in data: {length(common_genes)} (out of {nrow(gene_set)})", .sep = "\n"))
        
        if (length(common_genes) == 0) { 
            next
        }
        
        ES.vector <- vector(length = ncol(df))
        
        for (sample in seq_along(1:ncol(df))) {
            
            # Arrange genes by INCREASING expression (remember, ranked previously)
            gene_list <- order(df[, sample], decreasing = TRUE)   
            ordered_genes <- df[gene_list, sample]
            
            hit <- names(ordered_genes) %in% common_genes
            no_hit <- 1 - hit 
            
            ordered_genes <- abs(ordered_genes)^0.25
            
            sum.correl <- sum(ordered_genes[hit == 1])
            
            F0 <- 
                (no_hit/sum(no_hit)) |> 
                cumsum()
            
            Fn <- 
                ((hit*ordered_genes)/sum.correl) |> 
                cumsum()
            
            ES <- (Fn - F0) |> 
                sum()
            
            ES.vector[sample] <- ES
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