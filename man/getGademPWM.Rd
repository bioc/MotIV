 
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
path <- system.file(package="MotIV")
jaspar <- readPWMfile(paste(path,"/extdata/jaspar2010.txt",sep=""))
jaspar.scores <- readDBScores(paste(path,"/extdata/jaspar2010_PCC_SWU.scores",sep=""))

#####Input#####
data(FOXA1_rGADEM)
motifs <- getGademPWM(gadem)
motifs.trimed <- trimPWMedge(motifs, threshold=1)
}
\keyword{misc}
