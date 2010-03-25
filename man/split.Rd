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
path <- system.file("extdata/",package="MotIV")
jaspar <- readPWMfile(paste(path,"jasparDB2009.txt",sep=""))
jaspar.scores <- readDBScores(paste(path,"jaspar2009_PCC_SWU.scores",sep=""))

#####Input#####
load(paste(system.file("data/",package="MotIV"),"gadem_STAT1_MACS.rda",sep=""))
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

