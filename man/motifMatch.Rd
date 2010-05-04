\name{motifMatch}
\alias{motifMatch}
\title{Motifs Matches Analysis} 

\description{
Search for motifs matches corresponding to PWM. 
} 
\usage{
  motifMatch(inputPWM, database=jaspar, DBscores=jaspar.scores, cc="PCC", align="SWU", top=5, go=1, ge=0.5)
}
\arguments{
	\item{inputPWM}{A list of PWM.}
	\item{database}{A list of PWM corresponding to the database.}
	\item{DBscores}{A matrix object containing the scores associated to the database.}
	\item{cc}{The metric name to be used}
	\item{align}{The Alignment method to be used.}
	\item{top}{The number of identified transcription factors per motif.}
	\item{go}{Gap open penality.}
	\item{ge}{Gap extension penality.}
}
\details{
	For a set of PWMs given by \code{inputPWM}, this function realizes alignments with each motif of the \code{database} and returns the \code{top} best motifs.
	If no \code{database} is provided, the function will use jaspar by loading data(jaspar2010).
	If no \code{DBscores} is given, jaspar.scores from data(jaspar2010_scores) will be used.
	
	The e-value is computed according the matric name \code{cc} and is corrected by the \code{DBscores}. 
}
\value{A \code{motiv} object.}
\author{Eloi Mercier <\email{eloi.mercier@ircm.qc.ca}>}
\references{S Mahony, PE Auron, PV Benos, \code{DNA familial binding profiles made easy: comparison of various motif alignment and clustering strategies}, \emph{PLoS Computational Biology }(2007) 3(3):e61  }
\seealso{\code{generateDBScores}}
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
summary(foxa1.analysis.jaspar )
}
\keyword{misc}
