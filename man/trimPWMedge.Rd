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

\author{Eloi Mercier <\email{emercier@chibi.ubc.ca}>}

\seealso{
	 \code{makePWM}
}
\examples{
#####Database and Scores#####
path <- system.file(package="MotIV")
jaspar <- readPWMfile(paste(path,"/extdata/jaspar2010.txt",sep=""))
jaspar.scores <- readDBScores(paste(path,"/extdata/jaspar2010_PCC_SWU.scores",sep=""))

#####Input#####
data(FOXA1_rGADEM)
motifs <- getPWM(gadem)
motifs.trimed <- lapply(motifs,trimPWMedge, threshold=1)
}
\keyword{misc}
