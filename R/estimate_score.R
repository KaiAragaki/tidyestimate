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
#'
#' @examples

estimate_score <- function(df, is_affymetrix, tidy) {

    if(tidy) {
        df <- as.data.frame(df)
        rownames(df) <- df[, 1]
        df <- df[, -1]
    }
    
    # Sample rank normalization
    
    df <- df %>% 
        data.matrix() %>% 
        apply(2, rank, ties.method = "average")
    
    df <- 10000*df/nrow(df)
    
    gs.names <- names(SI_geneset) 
    
    score.matrix <- matrix(NA_real_, nrow = 2, ncol = ncol(df))
    
    for (gs.i in 1:2) {
        gene.set <- gs[gs.i,]
        gene.overlap <- intersect(gene.set, gene.names)
        print(paste0(gs.i, "Gene set: ", gs.names[gs.i], 
                    "\nNumber in provided dataset: ", length(gene.overlap), " (", length(gene.overlap), "/", length(gene.set), ")"))
        
        if (length(gene.overlap) == 0) { 
            next
        } else {
            ES.vector <- vector(length = Ns)
            
            # Enrichment score
            for (S.index in 1:Ns) {
                
                # Rank genes by expression
                gene.list <- order(m[, S.index], decreasing = TRUE)   
                ordered_genes <- m[gene.list, S.index]
                
                # Return the position of the overlapping genes (that is, the gene
                # list minus the ones that weren't there) in the supplied
                # dataframe
                gene.set2 <- match(gene.overlap, gene.names)
                
                TAG <- gene.list %in% gene.set2   # 1 (TAG) & 0 (no.TAG)
                # Inverse vector - 1 is 0, 0 is 1.
                no.TAG <- 1 - TAG 

                # Number of non-signature genes
                Nm <- length(gene.list) - length(gene.set2)
                
                # Unsure why there's an abs here - should all be positive
                # Too scared to remove it
                ordered_genes <- abs(ordered_genes)^0.25
                
                # Sum of the ranks for which there were hits
                sum.correl <- sum(ordered_genes[TAG == 1])
                
                # Non-hits divided by number of non-signature genes
                # A penalty is therefore 1/number of non-signature genes
                P0 <- no.TAG/Nm
                
                # Running sum of non-hits
                # c(0.01, 0.01, 0.02, 0.03, 0.04, 0.04, 0.05...)
                F0 <- cumsum(P0)
                
                # The proportion each gene contributes (sum should equal 1)
                # c(0, 0.3, 0, 0, 0, 0.12, 0...)
                Pn <- TAG * ordered_genes / sum.correl
                
                # c(0, 0.3, 0.3, 0.3, 0.3, 0.42, 0.42...)
                Fn <- cumsum(Pn)
                
                # A running enrichment score
                # c(-0.01, 0.29, 0.28, 0.27, 0.26, 0.38, 0.37...)
                RES <- Fn - F0
                ES.vector[S.index] <- sum(RES)
            }
            score.matrix[gs.i, ] <- ES.vector
        }
    }

    # PICK UP HERE
    
    score.data <- data.frame(score.matrix)
    names(score.data) <- sample.names
    row.names(score.data) <- gs.names
    estimate.score <- apply(score.data, 2, sum)
   
    if (!is_affymetrix){
        score.data <- rbind(score.data, estimate.score)
        rownames(score.data) <- c("StromalScore",
                                  "ImmuneScore",
                                  "ESTIMATEScore")
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
            if (est_i >= 0) {
                next
            } else {
                message(paste(sample.names[i],": out of bounds", sep=""))
            }
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
    outputGCT(score.data, output.ds)
}

# The ideal output format would have four columns. One would be sample name, the
# other three would be Stromal, Immune, and ESTIMATE score. I guess Affy gets a
# tumor purity score as well.