% $Id: e00-estimate-package.Rd 10 2016-09-28 18:15:01Z proebuck $

\name{estimate-package}
\alias{estimate-package}
\alias{estimate}
\docType{package}
\title{ESTIMATE}
\description{
  ESTIMATE (Estimation of STromal and Immune cells in MAlignant Tumor tissues
  using Expression data) is a package for predicting tumor purity, and the
  presence of stromal/immune cells in tumor tissues using gene expression data.
}
\details{
  \tabular{ll}{
    Package: \tab estimate\cr
    Type:    \tab Package\cr
    Version: \tab \Sexpr[results=text]{packageDescription("estimate")$Version}\cr
    Date:    \tab \Sexpr[results=text]{packageDescription("estimate")$Date}\cr
    License: \tab \Sexpr[results=text]{packageDescription("estimate")$License}\cr
  }
  For a complete list of functions, use \code{library(help="estimate")}.
}
\references{
  Subramanian A, et al.\cr
  "\href{http://www.pnas.org/content/102/43/15545.abstract}{Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles.}"\cr
  Proc Natl Acad Sci U S A 2005, 102:15545-15550.\cr
\cr
  Barbie DA, et al.\cr
  "\href{http://www.nature.com/nature/journal/v462/n7269/abs/nature08460.html}{Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1.}"\cr
  Nature 2009, 462:108-112.  \cr
\cr
  Verhaak RG, et al.\cr
  "\href{http://www.cell.com/cancer-cell/abstract/S1535-6108(09)00432-2}{Integrated genomic analysis identifies clinically relevant subtypes of glioblastoma characterized by abnormalities in PDGFRA, IDH1, EGFR, and NF1.}"\cr
  Cancer Cell 2010, 17:98-110.  \cr
}
\author{
  Kosuke Yoshihara \email{kyoshihara@mdanderson.org},
  P. Roebuck \email{proebuck@mdanderson.org}
}
\keyword{package}

