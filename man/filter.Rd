\name{filter}
\docType{methods}
\alias{filter}
\alias{filter-methods}
\alias{filter,motiv,list-method}
\alias{filter,motiv,filters-method}
\title{Filter Motifs } 

\description{
This function selects motifs according to a set of filters. 
} 
\usage{
\S4method{filter}{motiv,filters}(x, f, exact=FALSE, verbose=TRUE)
}
\arguments{
	\item{x}{An object of class \code{motiv}.}
	\item{f}{A filter or a set of filter for \code{motiv} object.}
	\item{verbose}{If \code{FALSE}, no output will be print.}
	\item{exact}{If \code{TRUE}, search only for perfect name match.}
}
\details{
	This function is used to select motifs that correspond to the filters.
	
	Many \code{filter} could be pass in argument separated by coma. They will be considered independently.
}
\value{A \code{motiv} object.}
\author{Eloi Mercier <\email{emercier@chibi.ubc.ca}>}

\seealso{
	 \code{setFilter} , \code{split}, \code{combine}
}
\examples{
#####Database and Scores#####
path <- system.file(package="MotIV")
jaspar <- readPWMfile(paste(path,"/extdata/jaspar2010.txt",sep=""))
jaspar.scores <- readDBScores(paste(path,"/extdata/jaspar2010_PCC_SWU.scores",sep=""))

#####Input#####
data(FOXA1_rGADEM)
motifs <- getPWM(gadem)
motifs.trimed <- lapply(motifs,trimPWMedge, threshold=1)

#####Analysis#####
foxa1.analysis.jaspar <- motifMatch(inputPWM=motifs,align="SWU",cc="PCC",database=jaspar,DBscores=jaspar.scores,top=5)
summary(foxa1.analysis.jaspar )

#####Filters#####
f.foxa1<-setFilter(name="", tfname="FOXA1", top=3, evalueMax=10^-5)
f.ap1 <- setFilter (tfname="AP1", top=3)
f.foxa1.ap1 <- f.foxa1 | f.ap1 
foxa1.filter <- filter(foxa1.analysis.jaspar, f.foxa1.ap1, exact=FALSE, verbose=TRUE)
foxa1.split <- split(foxa1.analysis.jaspar, c(f.foxa1, f.ap1) , drop=FALSE, exact=FALSE, verbose=TRUE)
foxa1.filter.combine <- combineMotifs(foxa1.filter, c(f.foxa1, f.ap1), exact=FALSE, name=c("FOXA1", "AP1"), verbose=TRUE)
}
\keyword{methods}
