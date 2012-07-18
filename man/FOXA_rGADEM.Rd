\name{FOXA1_rGADEM}
\docType{data}
\alias{FOXA1_rGADEM}
\alias{FOXA}
\alias{gadem}
\title{Dataset for FOXA1} 

\description{
This dataset contains results obtained by rGADEM for the FOXA1 transcription factor.
} 
\usage{
gadem
}

\references{\url{http://genomebiology.com/2008/9/9/R137}}
\examples{
#####Database and Scores#####
path <- system.file(package="MotIV")
data(jaspar2010)
data(jaspar2010_scores)

#####Input#####
data(FOXA1_rGADEM)
motifs <- getPWM(gadem)
motifs.trimed <- lapply(motifs,trimPWMedge, threshold=1)

#####Analysis#####
foxa1.analysis.jaspar <- motifMatch(inputPWM=motifs,align="SWU",cc="PCC",database=jaspar,DBscores=jaspar.scores,top=5)
summary(foxa1.analysis.jaspar )
}
\keyword{misc}
