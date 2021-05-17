% $Id: e90-SI_geneset-data.Rd 10 2016-09-28 18:15:01Z proebuck $

\name{SI_geneset-data}
\alias{SI_geneset}
\docType{data}
\title{two signatures for estimate}
\description{
  This data set contains two gene signatures (stromal and immune signatures).
  The stromal signature including 141 stroma-specific genes is designed to
  capture the presence of stroma in tumor tissue. The immune signature
  consisting of 141 immune cell-specific genes represents the infiltration of
  immune cells in tumor tissue.
}
\usage{
data(SI_geneset)
}
\format{
  The object \code{SI_geneset} is a \code{data.frame} with components:
  \describe{
    \item{StromalSignature}{character string specifying 141 genes in stromal
      signature}
    \item{ImmuneSignature}{character string specifying 141 genes in immune
      signature}
  }
}
\note{
  This data set was based on \file{SI_geneset.gmt} (in the \file{extdata}
  directory).
}
\author{
  Kosuke Yoshihara \email{kyoshihara@mdanderson.org}
}
\keyword{dataset}

