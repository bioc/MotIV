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
path <- system.file("extdata/",package="MotIV")
jaspar <- readPWMfile(paste(path,"jasparDB2009.txt",sep=""))
jaspar.scores <- readDBScores(paste(path,"jaspar2009_PCC_SWU.scores",sep=""))

#####Input#####
load(paste(system.file("data/",package="MotIV"),"gadem_STAT1_MACS.rda",sep=""))
motifs <- getGademPWM(gadem)
motifs.trimed <- trimPWMedge(motifs, threshold=1)

#####Analysis#####
stat1.analysis.jaspar <- motifMatch(inputPWM=motifs,align="SWU",cc="PCC",database=jaspar,DBscores=jaspar.scores,top=5)
viewMotifs(stat1.analysis.jaspar, 5)
}
\keyword{misc}
