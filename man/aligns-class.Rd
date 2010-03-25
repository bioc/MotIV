\name{alignments-class}
\Rdversion{1.1}
\docType{class}
\alias{alignments-class}

\title{Class "alignments" }
\description{This object contains the alignments found by a MotIV analysis.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("alignments", TF, evalue, sequence, match, strand )}.
}
\section{Slots}{
\describe{
	 \item{TF}{Object of class \code{"TF"}  }
    \item{evalue}{The e-value of the alignment. }
    \item{sequence}{The input sequence aligned.}
    \item{match}{The TF sequence wich as been matched. }
    \item{strand}{The strand of the alignment.}
}}

\author{Eloi Mercier <\email{eloi.mercier@ircm.qc.ca}>}
\seealso{
	 \code{motiv} , \code{matches}, \code{transcriptionFactor}
}
\keyword{classes}