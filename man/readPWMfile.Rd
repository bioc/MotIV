 
 
\name{readPWMfile}
\alias{readPWMfile}
\title{Read Transfac File} 

\description{
This function is use to read standard Transfac type file.} 
\usage{
  readPWMfile(file)
}    
\arguments{
	\item{file}{Transfac file's name.}
	}
\details{This function is designed to read standard Transfac type file. For more information about the format, please refere to \url{http://mcast.sdsc.edu/doc/transfac-format.html}}
\author{Eloi Mercier <\email{eloi.mercier@ircm.qc.ca}>}
\value{A list of matrix.}
\examples{
#####Database and Scores#####
path <- system.file(package="MotIV")
jaspar <- readPWMfile(paste(path,"/extdata/jaspar2010.txt",sep=""))
jaspar.scores <- readDBScores(paste(path,"/extdata/jaspar2010_PCC_SWU.scores",sep=""))
}
\keyword{misc}
 
