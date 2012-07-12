\name{as.data.frame}
\docType{methods}
\alias{as.data.frame}
\alias{as.data.frame-methods}
\alias{as.data.frame,motiv-method}
\title{Coerce to a Data Frame} 

\description{
This function coerces a MotIV object into a data frame. 
} 
\usage{
\S4method{as.data.frame}{motiv}(x)
}
\arguments{
	\item{x}{An object of class \code{motiv}.}
}
\details{
	'as.data.frame' returns a data frame. 
	
	This object regroups all the TF identified by \code{MotIV} with the corresponding evalue and alignments.
}
\value{A \code{data.frame} object.}
\author{Eloi Mercier <\email{emercier@chibi.ubc.ca}>}

\seealso{
	 \code{data.frame} , \code{viewAlignments}
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
foxa1.analysis.jaspar.df = as.data.frame(foxa1.analysis.jaspar)
head(foxa1.analysis.jaspar.df)
}
\keyword{methods}
