\name{viewMotifs-methods}
\docType{methods}
\alias{viewMotifs}
\alias{viewMotifs-methods}
\alias{viewMotifs,motiv-method}

\title{Print Identified Motifs} 

\description{
This function return a list of the identified motifs contained in a \code{motiv} object. 
} 
\usage{
	\S4method{viewMotifs}{motiv}(x,n=100)
}
\arguments{
	\item{x}{An object of class motiv (usualy provied by motifMatch).}
	\item{n}{The number of motifs shown.}
}
\details{This function shows the number of identified motif.}
\value{A list of motifs names.}

\author{Eloi Mercier <\email{eloi.mercier@ircm.qc.ca}>}
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
viewMotifs(foxa1.analysis.jaspar, 5)
}
\keyword{misc}
