% $Id: e90-common_genes-data.Rd 10 2016-09-28 18:15:01Z proebuck $

\name{common_genes-data}
\alias{common_genes}
\docType{data}
\title{10,412 common genes}
\description{
  This data set contains the information on 10,412 genes that are common among
  six platforms:
  \itemize{
    \item Affymetrix HG-U133Plus2.0
    \item Affymetrix HT-HG-U133A
    \item Affymetrix Human X3P
    \item Agilent 4x44K (G4112F)
    \item Agilent G4502A
    \item Illumina HiSeq RNA sequence
  }
  All gene information (according to NCBI) was current as of this package's
  release date.
}
\usage{
data(common_genes)
}
\format{
  The object \code{common_genes} is a \code{data.frame} with components:
  \describe{
    \item{EntrezID}{integer scalar specifying Entrez Gene ID}
    \item{GeneSymbol}{character string specifying official HUGO gene symbol}
    \item{Synonyms}{character string specifying alternative gene symbols,
      each delimited by vertical bar}
    \item{GeneName}{character string specifying official HUGO gene name}
    \item{Chromosome}{character string specifying chromosomal location}
  }
}
\author{
  Kosuke Yoshihara \email{kyoshihara@mdanderson.org}
}
\keyword{dataset}

