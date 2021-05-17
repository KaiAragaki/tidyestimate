% $Id: e50-outputGCT.Rd 10 2016-09-28 18:15:01Z proebuck $

\name{outputGCT}
\alias{outputGCT}
\title{Write gene expression data in GCT format}
\description{
  This function writes an input file in GCT format.
}
\usage{
outputGCT(input.f, output.f)
}
\arguments{
  \item{input.f}{character string specifying name of input file, or data frame
    including gene-level expression data}
  \item{output.f}{character string specifying name of output file}
}
\details{
  This function writes gene-level expression data to GCT format file. Input
  file must be a text format file.
}
\references{
  Reich M, et al.\cr
  "\href{http://www.nature.com/ng/journal/v38/n5/full/ng0506-500.html}{
    GenePattern 2.0.}"\cr
  Nature Genetics 2006, 38:500-501.\cr
}
\author{
  Kosuke Yoshihara \email{kyoshihara@mdanderson.org}
}
\examples{
in.file <- system.file("extdata", "sample_input.txt", package="estimate")
out.file <- tempfile(pattern="estimate", fileext=".gct")
outputGCT(in.file, out.file)
}
\keyword{convert}

