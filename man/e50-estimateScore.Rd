% $Id: e50-estimateScore.Rd 10 2016-09-28 18:15:01Z proebuck $

\name{estimateScore}
\alias{estimateScore}
\title{Calculation of stromal, immune, and ESTIMATE scores}
\description{
  This function computes stromal, immune, and ESTIMATE scores per sample using
  gene expression data.
}
\usage{
estimateScore(input.ds,
              output.ds,
              platform = c("affymetrix", "agilent", "illumina"))
}
\arguments{
  \item{input.ds}{character string specifying name of input GCT file containing
    stromal, immune, and estimate scores for each sample}
  \item{output.ds}{character string specifying name of output file}
  \item{platform}{character string indicating platform type. Defaults to
    "affymetrix"}
}
\details{
  This method is based on single sample gene set enrichment analysis (ssGSEA)
  algorithm. This function computes stromal, immune, and ESTIMATE scores using
  gene-level expression data. For Affymetrix platform data, tumor purity are
  derived from ESTIMATE scores by applying non-linear squares methods to
  \acronym{TCGA} Affymetrix expression data (n=995).
}
\value{
  Returns \code{data.frame} with components:
  \itemize{
    \item{StromalScore}{numeric scalar specifying the presence of stromal cells
      in tumor tissue}
    \item{ImmuneScore}{numeric scalar specifying the level of infiltrating
      immune cells in tumor tissue}
    \item{ESTIMATEScore}{numeric scalar specifying tumor cellularity}
    \item{TumorPurity}{numeric scalar specifying ESTIMATE-based tumor purity
      with value in range[0,1]}
  }
}
\references{
  Subramanian A, et al.\cr
  "\href{http://www.pnas.org/content/102/43/15545.abstract}{Gene set enrichment
    analysis: a knowledge-based approach for interpreting genome-wide expression
    profiles.}"\cr
  Proc Natl Acad Sci U S A 2005, 102:15545-15550.\cr
\cr
  Barbie DA, et al.\cr
  "\href{http://www.nature.com/nature/journal/v462/n7269/abs/nature08460.html}{
    Systematic RNA interference reveals that oncogenic KRAS-driven cancers
    require TBK1.}"\cr
  Nature 2009, 462:108-112.  \cr
\cr
  Verhaak RG, et al.\cr
  "\href{http://www.cell.com/cancer-cell/abstract/S1535-6108(09)00432-2}{
    Integrated genomic analysis identifies clinically relevant subtypes of
    glioblastoma characterized by abnormalities in PDGFRA, IDH1, EGFR, and
    NF1.}"\cr
  Cancer Cell 2010, 17:98-110.\cr
\cr
  Carter SL, et al.\cr
  "\href{http://www.nature.com/nbt/journal/v30/n5/full/nbt.2203.html?WT.ec_id=NBT-201205}{Absolute quantification of somatic DNA alterations in human cancer.}"\cr
  Nat Biotechnol 2012, 30:413-421.\cr
}
\author{
  Kosuke Yoshihara \email{kyoshihara@mdanderson.org}
}
\examples{
in.file <- system.file("extdata", "sample_input.gct", package="estimate")
out.file <- tempfile(pattern="estimate", fileext=".gct")
estimateScore(in.file, out.file)
}
\keyword{predict}

