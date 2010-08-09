\name{occurences}
\alias{cooccurences}
\alias{occurences}
\title{Motifs Occurences and Co-occurences} 

\description{
Get the number of motifs occurences and co-occurences from a \code{rGADEM} object.
} 
\usage{
occurences(gadem)
cooccurences(x)
}
\arguments{
	\item{gadem}{An object of type \code{rGADEM}.}
	\item{x}{A contingency table.}
}
\value{\code{occurences} returns the contingency table of the number of motifs per sequences.

This object can be put in \code{cooccurences} to return the number of sequences where two motifs appear together.
}

\author{Eloi Mercier <\email{eloi.mercier@ircm.qc.ca}>}

\examples{
data("FOXA1_rGADEM")
oc <- occurences (gadem)
coc <- cooccurences (oc)
coc
}
\keyword{misc}
