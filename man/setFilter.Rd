\name{setFilter}
\alias{setFilter}
\title{Set Motif Filter} 

\description{
This function is use to set a motif filter. 
} 
\usage{
  setFilter(name="", tfname="",evalueMax=1, top=10, lengthMax=100,valid=NULL)
}
\arguments{
	\item{name}{A name or a list of names.}
	\item{tfname}{A transcription factor name or a list of TF names.}
	\item{evalueMax}{An evalue between 0 and 1.}
	\item{top}{Defines the depth of the filter.}
	\item{lengthMax}{The maximum motif length.}
	\item{valid}{The alignment that should be considered as valid.}
}
\value{A \code{filter} object.}

\author{Eloi Mercier <\email{eloi.mercier@ircm.qc.ca}>}

\seealso{
	 \code{filter} , \code{split}, \code{combine}
}
\examples{
#####Database and Scores#####
path <- system.file(package="MotIV")
jaspar <- readPWMfile(paste(path,"/extdata/jaspar2010.txt",sep=""))
jaspar.scores <- readDBScores(paste(path,"/extdata/jaspar2010_PCC_SWU.scores",sep=""))

#####Input#####
data(FOXA1_rGADEM)
motifs <- getGademPWM(gadem)
motifs.trimed <- trimPWMedge(motifs, threshold=1)

#####Analysis#####
foxa1.analysis.jaspar <- motifMatch(inputPWM=motifs,align="SWU",cc="PCC",database=jaspar,DBscores=jaspar.scores,top=5)
summary(foxa1.analysis.jaspar )

#####Filters#####
f.foxa1<-setFilter(name="", tfname="FOXA1", top=3, evalueMax=10^-5)
f.ap1 <- setFilter (tfname="AP1", top=3)
f.foxa1.ap1 <- f.foxa1 | f.ap1 
foxa1.filter <- filter(foxa1.analysis.jaspar, f.foxa1.ap1, exact=FALSE, verbose=TRUE)
foxa1.split <- split(foxa1.analysis.jaspar, c(f.foxa1, f.ap1) , drop=FALSE, exact=FALSE, verbose=TRUE)
foxa1.filter.combine <- combine(foxa1.filter, c(f.foxa1, f.ap1), exact=FALSE, name=c("FOXA1", "AP1"), verbose=TRUE)
}
\keyword{misc}
