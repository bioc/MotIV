\name{filter-class}
\Rdversion{1.1}
\docType{class}
\alias{filter-class}


\title{Class "filter"}
\description{This object information to be apply as filter.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{ new("filter",name,tfname, top,evalueMax, lengthMax, valid)}.
}

\section{Slots}{
\describe{
	\item{name}{A name or a list of names.}
	\item{tfname}{A transcription factor name or a list of TF names.}
	\item{evalueMax}{An e-value between 0 and 1.}
	\item{top}{Defined the depth of the filter.}
	\item{lengthMax}{The maximum motif length.}
	\item{valid}{The alignment that should be considered as valid.}
	}
}

\details{
		This class \code{filter} is used to selected \code{motiv} objects according the filter's arguments.
}
\author{Eloi Mercier <\email{emercier@chibi.ubc.ca}>}

\seealso{
	 \code{setFilter} , \code{filter}, \code{split}, \code{combine}
}
\examples{
showClass("filter")
}
\keyword{classes}