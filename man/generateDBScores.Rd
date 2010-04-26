\name{generateDBScores }
\alias{generateDBScores}
\alias{readDBScores}
\alias{writeDBScores}
\title{Database Scores Functions} 

\description{
This functions are used to generate scores of a PWM database.
} 
\usage{
generateDBScores(inputDB,cc="PCC",align="SWU",nRand=1000,go=1,ge=0.5)
readDBScores(file)
writeDBScores(x, file)
}
\arguments{
	\item{inputDB}{A list of PWM corresponding to the database.}
	\item{cc}{The metric name to be used :}
	\item{align}{The Alignment method to be used.}
	\item{go}{Gap open penality.}
	\item{ge}{Gap extension penality.}
	\item{nRand}{The number of random PWM to be generated. The more higer it is, the more accurate score will be.}
	\item{file}{A character string naming a file.}
	\item{x}{A numeric matrix corresponding to a score.}
}
\details{
The score reflects the biais of the database. It is used to compute more precisely e-value alignments.

		\code{generateDBScores} : Based on database properties (suchs as length, zero rate, invariant colums ), nRand matrix are generated. A score is calculated for each matrix length with the specified alignment method and metric.
		
		The \code{score} is associated to a database and a alignment method and metric so you don't have to generate it each time you use the same database. Use the \code{writeDBScores} and \code{readDBScores} instead.
	\code{readDBScores} : Read a score file.
	\code{writeDBScores} : Write a score file.
}
\value{A numeric matrix. Columns correspond respectively to the first matrix length, second matrix length, variance, mean, matrix number, distance min and max.}
\section{Warning }{ Because of each matrix is compare to each other, computing time is exponential. You should be aware of this fact before provided a high nRand. 5000 is a good time/accuracy rate choice.}
\author{Shaun Mahony, modified by Eloi Mercier <\email{eloi.mercier@ircm.qc.ca}>}
\seealso{'readDBScores', 'writeDBScores'}
\references{Sandelin,A. and Wasserman,W.W. (2004) \code{Constrained binding site diversity within families of transcription factors enhances pattern discovery bioinformatics}. J. Mol. Biol., 338, 207/215.
}
\examples{
#####Database and Scores#####
path <- system.file(package="MotIV")
jaspar <- readPWMfile(paste(path,"/extdata/jaspar2010.txt",sep=""))
#jaspar.scores <- generateDBScores(inputDB=jaspar,cc="PCC",align="SWU",nRand=1000)
#writeDBScores(jaspar.scores,paste(path,"/extdata/jaspar_PCC_SWU.scores",sep=""))
jaspar.scores <- readDBScores(paste(path,"/extdata/jaspar2010_PCC_SWU.scores",sep=""))
}
\keyword{misc}
 
