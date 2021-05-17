% $Id: e50-filterCommonGenes.Rd 10 2016-09-28 18:15:01Z proebuck $

\name{filterCommonGenes}
\alias{filterCommonGenes}
\title{Intersect input data with 10,412 common genes}
\description{
  This function unifies different number of genes per platform against 10,412
  common genes.
}
\usage{
filterCommonGenes(input.f, output.f, id=c("GeneSymbol", "EntrezID"))
}
\arguments{
  \item{input.f}{character string specifying name of input file, or data.frame
    including gene-level expression data}
  \item{output.f}{character string specifying name of output file}
  \item{id}{character string indicating which gene identifier to use when
    matching.}
}
\details{
  The number of genes in expression data is different for each platform and
  this difference influences the computational results of stromal and immune
  scores. To compare stromal, immune and ESTIMATE scores across different
  platforms or calculate ESTIMATE-based tumor purity using Affymetrix expression
  data, users need to unify the gene identifiers of the input data against the
  common genes. Afterwards, the results of the intersection are written to disk
  in GCT format.
}
\author{
  Kosuke Yoshihara \email{kyoshihara@mdanderson.org}
}
\seealso{
  \code{\link{outputGCT}}
}
\examples{
in.file <- system.file("extdata", "sample_input.txt", package="estimate")
out.file <- tempfile(pattern="estimate", fileext=".gct")
filterCommonGenes(in.file, out.file)
}
\keyword{convert}

