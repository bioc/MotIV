\name{jaspar2010}
\docType{data}
\alias{jaspar}
\alias{jaspar.scores}
\title{Jaspar 2010 Database} 

\description{
Jaspar database and Jaspar score.
} 
\usage{
jaspar
jaspar.scores
}
\details{
	Jaspar is a well-known transcription factor database. Version 2010 contents 130 non-redundant matrix of TF binding sites.
	
	The jaspar scores have been computed with Pearson Correlation Coefficient and Smith-Waterman Ungapped alignments.
}

\source{http://jaspar.genereg.net/}
\references{Albin Sandelin, Wynand Alkema, P\"ar Engstr\"om, Wyeth W. Wasserman and Boris Lenhard, \code{JASPAR: an open-access database for eukaryotic transcription factor binding profiles}, \emph{Nucleic Acids Research}(2003)}
\seealso{\code{generateDBscores}, \code{motifMatch}}
\examples{
#####Database and Scores#####
path <- system.file(package="MotIV")
data(jaspar2010)
data(jaspar2010_scores)

#####Input#####
data(FOXA1_rGADEM)
motifs <- getGademPWM(gadem)
motifs.trimed <- trimPWMedge(motifs, threshold=1)

#####Analysis#####
foxa1.analysis.jaspar <- motifMatch(inputPWM=motifs,align="SWU",cc="PCC",database=jaspar,DBscores=jaspar.scores,top=5)
summary(foxa1.analysis.jaspar )
}
\keyword{misc}
