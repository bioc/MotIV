\name{split-methods}
\docType{methods}
\alias{split-methods}
\alias{split,motiv,list-method}
\alias{split,motiv,filters-method}

\title{Split Motiv Object} 
\description{
This function splits a 'motiv' object according filters. 
} 

\usage{
	\S4method{split}{motiv,filters}(x, f, exact=FALSE,drop=FALSE,verbose=TRUE,...)
}
\arguments{
	\item{x}{An object of class motiv (usualy provied by motifMatch).}
	\item{f}{A filter or a set of filter for motiv object.}
	\item{drop}{If \code{TRUE}, no match motifs will be dropped.}
	\item{verbose}{If \code{FALSE}, no output will be print.}
	\item{exact}{If \code{TRUE}, search only for perfect name match.}
	\item{...}{Further potential arguments passed to methods.}
}
\details{
	This function is used to split motifs that correspond to the filters.
	
	Many \code{filter} could be pass in argument separated by coma. They will be considered independently (coma is considered as \code{OR}).
}
\value{A list of \code{motiv} object.}
\author{Eloi Mercier <\email{eloi.mercier@ircm.qc.ca}>}
\seealso{
	 \code{setFilter} , \code{filter}, \code{combine}
}
\examples{
#####Database and Scores#####
path <- system.file(package="MotIV")
jaspar <- readPWMfile(paste(path,"/extdata/jaspar2010.txt",sep=""))
jaspar.scores <- readDBScores(paste(path,"/extdata/jaspar2010_PCC_SWU.scores",sep=""))

#####Input#####
data(FOXA1_rGADEM)
motifs <- getPWM(gadem)
motifs.trimed <- trimPWMedge(motifs, threshold=1)

#####Analysis#####
foxa1.analysis.jaspar <- motifMatch(inputPWM=motifs,align="SWU",cc="PCC",database=jaspar,DBscores=jaspar.scores,top=5)
summary(foxa1.analysis.jaspar )

#####Filters#####
f.foxa1<-setFilter(name="", tfname="FOXA1", top=3, evalueMax=10^-5)
f.ap1 <- setFilter (tfname="AP1", top=3)
f.foxa1.ap1 <- f.foxa1 | f.ap1 
foxa1.filter <- filter(foxa1.analysis.jaspar, f.foxa1.ap1, exact=FALSE, verbose=TRUE)
foxa1.split <- split(foxa1.analysis.jaspar, c(f.foxa1, f.ap1) , drop=FALSE, exact=FALSE, verbose=TRUE)
foxa1.filter.combine <- combine(foxa1.filter, c(f.foxa1, f.ap1), exact=FALSE, name=c("FOXA1", "AP1"), verbose=TRUE)
}
\keyword{methods}

