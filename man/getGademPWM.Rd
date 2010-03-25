 
\name{getGademPWM}
\alias{getGademPWM}
\title{Recover PWM } 

\description{
This function selects the PWMs contained in an object of type \code{gadem}.} 
\usage{
  getGademPWM(y)
}
\arguments{
	\item{y}{A \code{gadem} object.}
	}
\value{A list of PWM.}
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
}
\keyword{misc}
