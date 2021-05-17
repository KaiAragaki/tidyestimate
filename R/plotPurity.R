###
### $Id: plotPurity.R 12 2016-09-28 18:25:54Z proebuck $
###


##-----------------------------------------------------------------------------
plotPurity <- function(scores,
                       samples="all_samples",
                       platform=c("affymetrix", "agilent", "illumina"),
                       output.dir="estimated_purity_plots") {

    ## Check arguments
    stopifnot((is.character(scores) && length(scores) == 1 && nzchar(scores)) ||
              (inherits(scores, "connection") && isOpen(scores, "r")))
    stopifnot(is.character(output.dir) && length(output.dir) == 1 && nzchar(output.dir))
    platform <- match.arg(platform)
    
    if (platform != "affymetrix"){
        stop("not implemented")  
    }

    ## Begin processing

    ##-------------------------------------------------------------------------
    get_estimates_df <- function(scores) {
     #estimate <- read.table(scores, skip=2, header=TRUE, row.names=1, sep="\t")
     estimate <- read.delim(scores, skip=2, row.names=1)
     as.data.frame(t(estimate[, -1]))
    }

    ##-------------------------------------------------------------------------
    convert_row_estimate_score_to_tumor_purity <- function(x) {
      stopifnot(is.numeric(x))
      cos(0.6049872018 + 0.0001467884*x)
    }
  
    ## Read ESTIMATE data file
    estimate.df <- get_estimates_df(scores)
    samplenames <- rownames(estimate.df)
    Affy.model <- PurityDataAffy
    pred.p <- Affy.model[, 5:7]
    est <- estimate.df[, 3]
    est.new <- estimate.df[, 4]

    ## Create output directory
    dir.create(output.dir)

    ## ESTIMATE based tumor purity in scatterplot with prediction interval
    message("Plotting tumor purity based on ESTIMATE score")

    max.af <- max(Affy.model$ESTIMATEScore)
    min.af <- min(Affy.model$ESTIMATEScore)
  
    if (samples[1] == "all_samples"){
         Num.S <- nrow(estimate.df)
    } else {
         Num.S <- as.numeric(length(samples))
    }
  
    for (i in 1:Num.S) {
         if(samples[1] =="all_samples"){
             samplename <- samplenames[i]
        } else {
             samplename <- samples[i]
        }
    
        png.filename <- file.path(output.dir, sprintf("%s.png", samplename))
        png(filename=png.filename, width=480, height=480)
 
        geMin <- est[i] >= min.af
        leMax <- est[i] <= max.af
        withinMinMax <- geMin && leMax

        xlim <- if (!withinMinMax) {
             ## Expands plot boundary
             adjustment <- 500    # Arbitrary
             if (geMin) {
                 from <- min.af
                 to   <- est[i] + adjustment
            } else {
                 from <- est[i] - adjustment
                 to   <- max.af
            }
            c(from, to)
        } else {
             NULL
        }

        plot(Affy.model$tumor.purity~Affy.model$ESTIMATEScore, Affy.model,
        main=samplename,
        type="n",
        xlab="ESTIMATE score",
        xlim=xlim,
        ylab="Tumor purity",
        ylim=c(0, 1))
        par(new=TRUE)
        points(Affy.model$ESTIMATEScore, Affy.model$tumor.purity, cex=0.75, col="lightgrey")
        if (withinMinMax) {
             ## Prediction interval
             matlines(Affy.model$ESTIMATEScore, pred.p, lty=c(1, 2, 2), col="darkgrey")
        } else {
             matlines(Affy.model$ESTIMATEScore, pred.p, lty=c(1, 2, 2), col="darkgrey")
             par(new=TRUE)
             curve(convert_row_estimate_score_to_tumor_purity,
             from, to, n=10000, col="grey", ylim=c(0, 1), xlab="", ylab="")
        }
        points(est[i], est.new[i], pch=19, cex=1.25)
        abline(h=est.new[i], col="black", lty=2)
        abline(v=est[i], col="black", lty=2)

        dev.off()
    }
}

