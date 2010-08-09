\name{filters-methods}
\docType{methods}
\alias{summary,filter-method}
\alias{summary,filters-method}
\alias{names,filter-method}
\alias{names,filters-method}
\alias{show,filter-method}
\alias{&,filters,filters-method}
\alias{|,filters,filters-method}

\title{Filters Methods}
\description{Methods for filters object}

\usage{
	 \S4method{summary}{filter}(object)
	 \S4method{summary}{filters}(object)
	 \S4method{names}{filter}(x)
	 \S4method{names}{filters}(x)	 
}

\arguments{
	\item{object}{An object of class \code{filter}.}
	\item{x}{An object of class \code{filter}.}

}

\author{Eloi Mercier <\email{eloi.mercier@ircm.qc.ca}>}

\seealso{
	 \code{setFilter} , \code{filter}, \code{split}, \code{combine}
}
\examples{
showClass("filter")
}
\keyword{classes}