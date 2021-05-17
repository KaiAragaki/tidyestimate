###
### $Id: zzz.R 10 2016-09-28 18:15:01Z proebuck $
###
### Hooks called as package is loaded, attached, detatched, and unloaded.
###


##
## Package/Namespace Hooks
##

##-----------------------------------------------------------------------------
.onAttach <- function(libname, pkgname) {
    verbose <- getOption("verbose")
    if (verbose) {
        local({
            libraryPkgName <- function(pkgname, sep="_") {
                unlist(strsplit(pkgname, sep, fixed=TRUE))[1]
            }
            packageDescription <- function(pkgname) {
                fieldnames <- c("Title", "Version")
                metafile <- file.path(libname, pkgname, "DESCRIPTION")
                meta <- as.list(read.dcf(metafile, fieldnames))
                names(meta) <- fieldnames
                return(meta)
            }

            meta <- packageDescription(pkgname)
            msg <- sprintf("%s, version %s",
                           meta$Title, meta$Version)
            packageStartupMessage(msg)
            msg <- sprintf("Type library(help=%s) to see package documentation",
                           libraryPkgName(pkgname))
            packageStartupMessage(msg)
        })
    }
}


##-----------------------------------------------------------------------------
.onLoad <- function(libname, pkgname) {
    datasetNames <- c(
        "common_genes",
        "PurityDataAffy",
        "SI_geneset"
    )

    ## Make R CMD check shutup about "no visible binding for global variable"
    globalVariables(datasetNames)

    ## Package uses the "LazyData" setting to load the datasets when needed
    invisible(NULL)
}

