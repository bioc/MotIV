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
\author{Eloi Mercier <\email{eloi.mercier@ircm.qc.ca}>}

\seealso{
	 \code{setFilter} , \code{split}, \code{combine}
}
\examples{
#####Database and Scores#####
path <- system.file(package="MotIV")
jaspar <- readPWMfile(paste(path,"/extdata/jasparDB2009.txt",sep=""))
jaspar.scores <- readDBScores(paste(path,"/extdata/jaspar2009_PCC_SWU.scores",sep=""))

#####Input#####
load(paste(path,"/data/gadem_STAT1_MACS.rda",sep=""))
motifs <- getGademPWM(gadem)
motifs.trimed <- trimPWMedge(motifs, threshold=1)

#####Analysis#####
stat1.analysis.jaspar <- motifMatch(inputPWM=motifs,align="SWU",cc="PCC",database=jaspar,DBscores=jaspar.scores,top=5)
summary(stat1.analysis.jaspar )

#####Filters#####
f.statx<-setFilter(name="", tfname="[Ss][Tt][Aa][Tt]", top=3, evalueMax=10^-5)
f.ap1 <- setFilter (tfname="AP1", top=3)
f.statx.ap1 <- f.statx | f.ap1 
stat1.filter <- filter(stat1.analysis.jaspar, f.statx.ap1, exact=FALSE, verbose=TRUE)
stat1.split <- split(stat1.analysis.jaspar, c(f.statx, f.ap1) , drop=FALSE, exact=FALSE, verbose=TRUE)
stat1.filter.combine <- combine(stat1.filter, c(f.statx, f.ap1), exact=FALSE, name=c("STATx", "AP1"), verbose=TRUE)
}
\keyword{methods}
