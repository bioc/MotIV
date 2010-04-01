 
 
\name{readGademPWMFile}
\alias{readGademPWMFile}
\title{Read Gadem File} 

\description{
This function is use to read a gadem file containing PWM.} 
\usage{
  readGademPWMFile(file)
}
\arguments{
	\item{file}{File's name.}
	}
\details{This function is made to read typicaly output file from Gadem (v1.2). Standard name is 'observedPWMs.txt'.}
\author{Eloi Mercier <\email{eloi.mercier@ircm.qc.ca}>}
\examples{
#####Database and Scores#####
path <- system.file(package="MotIV")
jaspar <- readPWMfile(paste(path,"/extdata/jasparDB2009.txt",sep=""))
jaspar.scores <- readDBScores(paste(path,"/extdata/jaspar2009_PCC_SWU.scores",sep=""))

#####Input#####
load(paste(path,"/data/gadem_STAT1_MACS.rda",sep=""))
motifs <- getGademPWM(gadem)
motifs.trimed <- trimPWMedge(motifs, threshold=1)
}
\keyword{misc}
