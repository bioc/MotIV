\name{matches-class}
\Rdversion{1.1}
\docType{class}
\alias{matches-class}

\title{Class "matches"}
\description{This object contains the name of the input motif and all the matches found.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("matches", name, aligns, similarity, valid)}.
}
\section{Slots}{
\describe{
    \item{name}{Motif name.}
	\item{aligns}{Alignments found by \code{motifMatch}.}
	\item{similarity}{The optional name given to the motif.}
	\item{valid}{The alignment that should be considered as valid.}
	}
}

\author{Eloi Mercier <\email{eloi.mercier@ircm.qc.ca}>}

\seealso{
	 \code{motiv} , \code{alignments}, \code{tf}
}
\examples{
showClass("matches")
}
\keyword{classes}
