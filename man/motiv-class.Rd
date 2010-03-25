\name{motiv-class}
\Rdversion{1.1}
\docType{class}
\alias{motiv-class}

\title{Class "motiv"}
\description{This object contains all informations about the motiv analysis.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("motiv",input, bestMatch,argv)}.
}

\section{Slots}{
\describe{
    \item{input}{List of input PWM.}
    \item{bestMatch}{Object of class \code{"matches"}.}
    \item{argv}{List of arguments used.}
	}
}
\author{Eloi Mercier <\email{eloi.mercier@ircm.qc.ca}>}

\seealso{
	 \code{matches} , \code{alignments}, \code{transcriptionFactor}
}
\examples{
showClass("motiv")
}
\keyword{classes}
