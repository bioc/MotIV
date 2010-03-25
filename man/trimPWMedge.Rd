\name{trimPWMedge}
\alias{trimPWMedge}
\title{Trim PWM edge} 

\description{
This function is use to cut edges with low information content.
} 
\usage{
  trimPWMedge(x, threshold=1)
}
\arguments{
	\item{x}{A matrix representing a PWM.}
	\item{threshold}{A transcription factor name or a list of TF names.}
}
\value{A \code{PWM}.}

\author{Eloi Mercier <\email{eloi.mercier@ircm.qc.ca}>}

\seealso{
	 \code{makePWM}
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
}
\keyword{misc}
