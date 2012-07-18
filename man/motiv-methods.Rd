\name{motiv-methods}
\Rdversion{1.1}
\docType{methods}
\alias{motiv}
\alias{summary,motiv-method}
\alias{summary,list-method}
\alias{show,motiv-method}
\alias{names,motiv-method}
\alias{names<-,motiv-method}
\alias{length,motiv-method}
\alias{[,motiv-method}
\alias{similarity}
\alias{similarity,motiv-method}
\alias{similarity,list-method}

\title{Motiv methods}
\description{Methods for motiv objects.}

\usage{
	 \S4method{summary}{motiv}(object) 
	 \S4method{names}{motiv}(x)
	 \S4method{length}{motiv}(x)
	 \S4method{similarity}{motiv}(x)
	 \S4method{[}{motiv}(x,i,j=ANY, bysim=TRUE, ...,exact=TRUE, ignore.case=FALSE, drop=FALSE)
}

\arguments{
	\item{object}{An object of class \code{motiv}.}
	\item{x}{An object of class \code{motiv}.}
	\item{i}{A string reprsenting a motif name.}
	\item{j}{NOT USED.}
	\item{bysim}{If \code{TRUE}, select by similarity name.}
	\item{...}{Further potential arguments passed to methods.}
	\item{ignore.case}{if FALSE, the pattern matching is case sensitive and if TRUE, case is ignored during matching}
	\item{exact}{If \code{TRUE}, search only for perfect name match.}
	\item{drop}{If \code{TRUE}, no match motifs will be dropped.}
}

\author{Eloi Mercier <\email{emercier@chibi.ubc.ca}>}

\seealso{
	 \code{matches} , \code{alignments}, \code{tf}
}
\examples{
showClass("motiv")
}
\keyword{classes}
