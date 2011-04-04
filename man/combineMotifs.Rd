\name{combineMotifs}
\docType{methods}
\alias{combineMotifs}
\alias{combineMotifs-method}
\alias{combineMotifs,motiv,list-method}
\alias{combineMotifs,motiv,filters-method}
\title{Combine Motifs } 

\description{
This function combines motifs according to a set of filters. 
} 
\usage{
\S4method{combineMotifs}{motiv,filters}(x, y, name=NULL,exact=TRUE,verbose=TRUE)
}
\arguments{
	\item{x}{An object of class \code{motiv}.}
	\item{y}{A filter or a set of filter.}
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
\author{Eloi Mercier <\email{emercier@chibi.ubc.ca}>}
\seealso{
	 \code{setFilter} , \code{filter}, \code{split}
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
foxa1.filter.combine <- combineMotifs(foxa1.filter, c(f.foxa1, f.ap1), exact=FALSE, name=c("FOXA1", "AP1"), verbose=TRUE)
}
\keyword{methods}
