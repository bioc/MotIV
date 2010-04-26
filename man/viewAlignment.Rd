 
\name{viewAlignments}
\alias{viewAlignments}
\title{Print Motifs Alignments} 

\description{
This function return a list of the alignments of a motiv object for each motif. 
} 
\usage{
  viewAlignments(x)
}
\arguments{
	\item{x}{An object of class motiv (usualy provied by motifMatch).}
}
\details{This function shows the alignements for each motif.}
\author{Eloi Mercier <\email{eloi.mercier@ircm.qc.ca}>}
\examples{
#####Database and Scores#####
path <- system.file(package="MotIV")
jaspar <- readPWMfile(paste(path,"/extdata/jaspar2010.txt",sep=""))
jaspar.scores <- readDBScores(paste(path,"/extdata/jaspar2010_PCC_SWU.scores",sep=""))

#####Input#####
data(FOXA1_rGADEM)
motifs <- getGademPWM(gadem)
motifs.trimed <- trimPWMedge(motifs, threshold=1)

#####Analysis#####
foxa1.analysis.jaspar <- motifMatch(inputPWM=motifs,align="SWU",cc="PCC",database=jaspar,DBscores=jaspar.scores,top=5)
summary(foxa1.analysis.jaspar )
viewAlignments(foxa1.analysis.jaspar )
}
\keyword{misc}
