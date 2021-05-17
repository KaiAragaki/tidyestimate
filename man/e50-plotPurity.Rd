% $Id: e50-plotPurity.Rd 10 2016-09-28 18:15:01Z proebuck $

\name{plotPurity}
\alias{plotPurity}
\title{Plot tumor purity}
\description{
  Plot tumor purity based on ESTIMATE score.
}
\usage{
plotPurity(scores,
           samples="all_samples",
           platform = c("affymetrix", "agilent", "illumina"),
           output.dir="estimated_purity_plots")
}
\arguments{
  \item{scores}{a GCT file output by \code{estimateScore} method containing
    stromal, immune, ESTIMATE scores, and tumor purity per each sample}
  \item{samples}{vector of character strings specifying sample names to be
    plotted. Defaults to "all_samples", which creates plots for all input
    samples}
  \item{platform}{character string indicating platform type. Defaults to
    "affymetrix"}
  \item{output.dir}{character string specifying a directory where plots should
    be saved. Defaults to "estimated_purity_plots"}
}
\details{
  This function produces scatterplots (as PNG files) for each requested sample;
  it plots tumor purity against ESTIMATE score. At present, only the Affymetrix
  platform is supported.
}
\author{
  Kosuke Yoshihara \email{kyoshihara@mdanderson.org}
}
\examples{
in.file <- system.file("extdata", "sample_purity.gct", package="estimate")
plotPurity(scores=in.file)
}
\keyword{figure}

