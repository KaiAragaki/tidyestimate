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
#' @export
#' @import 
#' @examples

estimate_score <- function(df, is_affymetrix, tidy) {

    
    # TODO:
    # Better naming
    # Reimplement sample name bit at the end
    # Reimpliment msg adding
    # Reimplement affy part
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

    gene_sets <- SI_geneset
    gene_set_names <- colnames(gene_sets)
    score.matrix <- matrix(NA_real_, nrow = 2, ncol = ncol(df))
    
    for (i in 1:2) {
        gene_set <- gene_sets[, i]
        sig_genes_in_df <- intersect(gene_set, rownames(df))
        
        if (length(sig_genes_in_df) == 0) { 
            next
        } else {
            ES.vector <- vector(length = ncol(df))
            
            # Enrichment score
            for (sample in seq_along(1:ncol(df))) {
                
                # Arrange genes by INCREASING expression (remember, ranked previously)
                gene_list <- order(df[, sample], decreasing = TRUE)   
                ordered_genes <- df[gene_list, sample]

                hit <- names(ordered_genes) %in% sig_genes_in_df

                # Inverse vector - 1 is 0, 0 is 1.
                no_hit <- 1 - hit 

                # Unsure why there's an abs here - should all be positive
                # Too scared to remove it
                ordered_genes <- abs(ordered_genes)^0.25
                
                # Sum of the ranks for which there were hits
                sum.correl <- sum(ordered_genes[hit == 1])
                
                # Running sum of non-hits.
                # Each non-hit incurs a penalty of 1/(# non-sig-genes)
                # c(0.01, 0.01, 0.02, 0.03, 0.04, 0.04, 0.05...)
                F0 <- 
                    (no_hit/sum(no_hit)) |> 
                    cumsum()

                # Running sum of the proportion each gene contributes 
                # (sum should equal 1)
                # c(0, 0.3, 0, 0, 0, 0.12, 0...)
                Fn <- 
                    ((hit*ordered_genes)/sum.correl) |> 
                    cumsum()

                # Sum of running ES
                # c(-0.01, 0.29, 0.28, 0.27, 0.26, 0.38, 0.37...)
                ES <- (Fn - F0) |> 
                    sum()
                ES.vector[sample] <- ES
            }
            score.matrix[i, ] <- ES.vector
        }
    }

    score.data <- data.frame(score.matrix)
    names(score.data) <- colnames(df)
    row.names(score.data) <- gene_set_names
    estimate.score <- apply(score.data, 2, sum)
   
    if (!is_affymetrix){
        score.data <- rbind(score.data, estimate.score)
        rownames(score.data) <- c("StromalScore",
                                  "ImmuneScore",
                                  "ESTIMATEScore")
        score.data
    } else {
        # Calculate ESTIMATE-based tumor purity (Affymetrix-specific)
        convert_row_estimate_score_to_tumor_purity <- function(x) {
            stopifnot(is.numeric(x))
            cos(0.6049872018 + 0.0001467884*x)
        }
        
        est.new <- NULL
        for (i in 1:length(estimate.score)) {
            est_i <- convert_row_estimate_score_to_tumor_purity(estimate.score[i])
            est.new <- rbind(est.new, est_i)
        }
        colnames(est.new) <- c("TumorPurity")
        estimate.t1 <- cbind(estimate.score, est.new)
        x.bad.tumor.purities <- estimate.t1[, "TumorPurity"] < 0
        estimate.t1[x.bad.tumor.purities, "TumorPurity"] <- NA
        score.data <- rbind(score.data, t(estimate.t1))
        rownames(score.data) <- c("StromalScore",
                                  "ImmuneScore",
                                  "ESTIMATEScore",
                                  "TumorPurity")
    }
}