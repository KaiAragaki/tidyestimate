###
### $Id: filterCommonGenes.R 21 2016-10-04 21:13:08Z proebuck $
###


##-----------------------------------------------------------------------------
filterCommonGenes <- function(input.f,
                              output.f,
                              id=c("GeneSymbol", "EntrezID")) {

    ## Check arguments
    stopifnot((is.character(input.f) && length(input.f) == 1 && nzchar(input.f)) ||
              (inherits(input.f, "connection") && isOpen(input.f, "r")))
    stopifnot((is.character(output.f) && length(output.f) == 1 && nzchar(output.f)))
    id <- match.arg(id)   
     
    ## Read input data
    input.df <- read.table(input.f,
                           header=TRUE,
                           row.names=1,
                           sep="\t", 
                           quote="",
                           stringsAsFactors=FALSE)
     
    merged.df <- merge(common_genes, input.df, by.x=id, by.y="row.names")
    rownames(merged.df) <- merged.df$GeneSymbol
    merged.df <- merged.df[, -1:-ncol(common_genes)]
    print(sprintf("Merged dataset includes %d genes (%d mismatched).",
                  nrow(merged.df),
                  nrow(common_genes) - nrow(merged.df)))
    outputGCT(merged.df, output.f)
}

