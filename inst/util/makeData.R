local({
    data.dir    <- system.file("data", package="estimate")
    extdata.dir <- system.file("extdata", package="estimate")

    ## 1. common_genes
    txtfile <- file.path(extdata.dir, "common_genes.txt")
    dataset <- file.path(data.dir, sprintf("%s.RData", basename(txtfile)))
    common_genes <- read.delim(txtfile,
                               header=TRUE,
                               quote="",
                               stringsAsFactors=FALSE)
    save(common_genes, file=dataset)

    ## 2. SI-geneset
    gmtfile <- file.path(extdata.dir, "SI_geneset.gmt")
    dataset <- file.path(data.dir, sprintf("%s.RData", basename(gmtfile)))
    SI_geneset <- read.delim(gmtfile,
                             header=FALSE,
                             row.names=1,
                             quote="",
                             stringsAsFactors=FALSE)
    save(SI_geneset, file=dataset)

    ## 3. PurityDataAffy
    txtfile <- file.path(extdata.dir, "PurityDataAffy.txt")
    dataset <- file.path(data.dir, sprintf("%s.RData", basename(txtfile)))
    PurityDataAffy <- read.delim(txtfile,
                                 header=TRUE,
                                 row.names=1,
                                 quote="",
                                 stringsAsFactors=FALSE)
    save(PurityDataAffy, file=dataset)
})

