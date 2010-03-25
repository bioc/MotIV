\name{combine}
\docType{methods}
\alias{combine}
\alias{combine-method}
\alias{combine,motiv,list-method}
\alias{combine,motiv,filters-method}
\title{Combine Motifs } 

\description{
This function combines motifs according to a set of filters. 
} 
\usage{
\S4method{combine}{motiv,filters}(x, f, name=NULL,exact=FALSE,verbose=TRUE)
}
\arguments{
	\item{x}{An object of class \code{motiv}.}
	\item{f}{A filter or a set of filter.}
	\item{name}{Name(s) to be given for similar motifs.}
	\item{verbose}{If \code{FALSE}, no output will be print.}
	\item{exact}{If \code{TRUE}, search only for perfect name match.}
}
\details{
	This function is used to consider some motifs as a unique motif or \code{similar} motifs.
	
	Many \code{filters} could be pass in argument separated by coma. They will be considered independently (coma is considered as \code{OR}).

	If a name or a vector of name is provided, it will be used to assign new name for similar motif  to the corresponding filter.
Else, a generic name is used. 
}
\value{A \code{motiv} object.}
\author{Eloi Mercier <\email{eloi.mercier@ircm.qc.ca}>}
\seealso{
	 \code{setFilter} , \code{filter}, \code{split}
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
stat1.filter.combine <- combine(stat1.filter, c(f.statx, f.ap1), exact=FALSE, name=c("STATx", "AP1"), verbose=TRUE)
}
\keyword{methods}
