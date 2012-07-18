
\name{motifDistances}
\alias{motifDistances}
\alias{motifHclust}
\alias{motifCutree}
\title{Clustering PWMs Computation} 

\description{
Set of functions to perfom clustering of PWMs.} 
\usage{
	motifDistances(inputPWM, DBscores=jaspar.scores, cc="PCC", align="SWU", top=5, go=1, ge=0.5) 
	motifHclust(x,...)
	motifCutree(tree,k=NULL, h=NULL)
}
\arguments{
	\item{inputPWM, DBscores, cc, align, top, go, ge}{Option for the PWMs distances computation. Refere to \code{motifMatch}.}
	\item{x,...}{Arguments to pass to the hclust function. See \code{hclust}.}
	\item{tree, k, h}{Arguments to pass to the cutree function. See \code{cutree}.}
}
\details{
This function are made to perform motifs clustering. 

The \sQuote{motifDistances} function computes the distances between each pair of motifs using the specified alignment.

The \sQuote{motifHclust} and \sQuote{motifCutree} functions are simple redefinition of \sQuote{hclust} and \sQuote{cutree}.
}
\author{Eloi Mercier <\email{emercier@chibi.ubc.ca}>}
\examples{
#####Database and Scores#####
path <- system.file(package="MotIV")
jaspar <- readPWMfile(paste(path,"/extdata/jaspar2010.txt",sep=""))
jaspar.scores <- readDBScores(paste(path,"/extdata/jaspar2010_PCC_SWU.scores",sep=""))

#####Input#####
data(FOXA1_rGADEM)
motifs <- getPWM(gadem)
motifs.trimed <- lapply(motifs,trimPWMedge, threshold=1)

#####Analysis#####
foxa1.analysis.jaspar <- motifMatch(inputPWM=motifs,align="SWU",cc="PCC",database=jaspar,DBscores=jaspar.scores,top=5)

#####Clustering#####
d <- motifDistances(getPWM(foxa1.analysis.jaspar))
hc <- motifHclust(d)
plot(hc)
f <- motifCutree(hc, k=2)
foxa1.combine <- combineMotifs(foxa1.analysis.jaspar, f, exact=FALSE, name=c("Group1", "Group2"), verbose=TRUE)

}

\keyword{misc}
 
