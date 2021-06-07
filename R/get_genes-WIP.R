get_aliases <- function(common_genes) {
        
        bm <- biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
        
        biomaRt::getBM(
                attributes = c('entrezgene_id', 'hgnc_symbol'),
                filters = 'entrezgene_id',
                values = common_genes$EntrezID, 
                mart = bm)
}
