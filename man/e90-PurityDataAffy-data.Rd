% $Id: e90-PurityDataAffy-data.Rd 10 2016-09-28 18:15:01Z proebuck $

\name{PurityDataAffy-data}
\alias{PurityDataAffy}
\docType{data}
\title{Affymetrix data}
\description{
  This data set contains stromal, immune, and ESTIMATE scores in all
  Affymetrix expression data (n=995) that was used to develop the formula for
  predicting tumor purity based on raw estimate score.

  This data set also contains tumor purity based on ABSOLUTE algorithm,
  predicted tumor purity, and 95% confidence intervals per sample based on
  ESTIMATE algorithm.
}
\usage{
data(PurityDataAffy)
}
\format{
  The object \code{PurityDataAffy} is a \code{data.frame} with components:
  \describe{
    \item{tumor.purity}{numeric scalar specifying tumor purity calculated by
      ABSOLUTE algorithm}
    \item{StromalScore}{numeric scalar specifying the presence of stromal cells
      in tumor tissue}
    \item{ImmuneScore}{numeric scalar specifying the level of infiltrating
      immune cells in tumor tissue}
    \item{ESTIMATEScore}{numeric scalar specifying tumor cellularity}
    \item{fit}{numeric scalar specifying estimated tumor purity based on
      ESTIMATE algorithm}
    \item{lwr.p}{numeric scalar specifying 5 percent confidence interval}
    \item{upr.p}{numeric scalar specifying 95 percent confidence interval}
  }
}
\source{
  The Cancer Genome Atlas
}
\references{
  Carter SL, et al.\cr
  "\href{http://www.nature.com/nbt/journal/v30/n5/full/nbt.2203.html?WT.ec_id=NBT-201205}{Absolute quantification of somatic DNA alterations in human cancer.}"\cr
  Nat Biotechnol 2012, 30:413-421.\cr
}
\author{
  Kosuke Yoshihara \email{kyoshihara@mdanderson.org}
}
\keyword{dataset}

