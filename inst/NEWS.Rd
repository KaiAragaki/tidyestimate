% $Id: NEWS.Rd 20 2016-10-04 00:34:04Z proebuck $
\name{NEWS}
\title{News for Package 'estimate'}
\section{Changes in estimate version 1.0.13}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Now available from R-Forge!
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Denoted package imports.
      \item Revamped handling of internal datasets, lazyloading them
        into the package namespace.
      \item Squashed multiple R CMD check warnings.
    }
  }
}
\section{Changes in estimate version 1.0.12}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Unknown - changes were not logged in version control.
    }
  }
}
\section{Changes in estimate version 1.0.11}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Add "util" folder into "inst" directory and
        	put makeData.R into this directory.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item I fixed the places of ESTIMATE_vignette.Rnw and
        	ESTIMATE_vignette.pdf.
      \item Rewrite of documentation properly.
    }
  }
}
\section{Changes in estimate version 1.0.10}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Replace unifyCommonGenes.R with
        	filterCommonGenes.R, which includes
        	new argument "id".
      \item Revise data(common_genes)
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Rewrite of documentation properly.
      \item Code cleanup in outputGCT.R and estimateScore.R.
    }
  }
}
\section{Changes in estimate version 1.0.9}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Unify estimateScore.R and estimateScoreyAffy.R to 
        	new "estimateScore.R" by adding new argument "platform".
      \item Replase purityPlotAffy.R with purityPlot.R by adding
        	new argument "platform".
      \item Replace GeneSymbols and GeneName with HUGO official
        	Gene Symbols and Gene names in data(common_genes).
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Rewrite of documentation properly.
      \item Code cleanup.
    }
  }
}
\section{Changes in estimate version 1.0.8}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Added NEWS.Rd.
      \item Added ESTIMATE_vignette.Rnw.
      \item Replace CalculatePurityAffy.R with estimateScoreAffy.R.
      \item Added current NCBI Official Symbol, Synonyms, Description, 
        	and chromosome per each gene to data(common_genes).
      \item Added new argument "Sample.N" into plotPurityAffy.R.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Corrected code responsible for possible problems; 
        	plotPurityAffy: no visible binding for global variable 'x'.
      \item Revised misdescription in documentation for plotPurityAffy.
      \item Rewrite of documentation properly.
      \item Code cleanup.
    }
  }
}
\section{Changes in estimate version 1.0.7}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Added NCBI Entrez Gene ID per each gene to data(common_genes).
    }
  }
}
\section{Changes in estimate version 1.0.6}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Revised calculatePurityAffy and plotPurityAffy to deal with
        	extremely high estimate score.
      \item Excluded "RColorBrewer", "gtools", "vertification", "survival"
        	packages from "Depends".
    }
  }
}
\section{Changes in estimate version 1.0.5}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Added additional information to data(PurityDataAffy).
      \item Added new convert function instead of old one.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Massive rewrite of documentation.
      \item Massive code cleanup.
    }
  }
}
\section{Changes in estimate version 1.0.4}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Modified estimateScore to be able to run without 
        	GSEA.EnrichmentScore5.R.
      \item Rename R-script file names to easily understand 
        	the function of each script.
    }
  }
}
\section{Changes in estimate version 1.0.3}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Deleted estimate.OPAM.R.
      \item Added new R-scripts instead of estimate.OPAM.R.
    }
  }
}
\section{Changes in estimate version 1.0.2}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Massive rewrite of documentation.
      \item Massive code cleanup.
      \item Actually passes R CMD check.
    }
  }
}
