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
path <- system.file("extdata/",package="MotIV")
jaspar <- readPWMfile(paste(path,"jasparDB2009.txt",sep=""))
jaspar.scores <- readDBScores(paste(path,"jaspar2009_PCC_SWU.scores",sep=""))

#####Input#####
load(paste(system.file("data/",package="MotIV"),"gadem_STAT1_MACS.rda",sep=""))
motifs <- getGademPWM(gadem)
motifs.trimed <- trimPWMedge(motifs, threshold=1)

#####Analysis#####
stat1.analysis.jaspar <- motifMatch(inputPWM=motifs,align="SWU",cc="PCC",database=jaspar,DBscores=jaspar.scores,top=5)
summary(stat1.analysis.jaspar )

#####Filters#####
f.statx<-setFilter(name="", tfname="[Ss][Tt][Aa][Tt]", top=3, evalueMax=10^-5)
f.ap1 <- setFilter (tfname="AP1", top=3)
f.statx.ap1 <- f.statx | f.ap1 
stat1.filter <- filter(stat1.analysis.jaspar, f.statx.ap1, exact=FALSE, verbose=TRUE)
stat1.split <- split(stat1.analysis.jaspar, c(f.statx, f.ap1) , drop=FALSE, exact=FALSE, verbose=TRUE)
stat1.filter.combine <- combine(stat1.filter, c(f.statx, f.ap1), exact=FALSE, name=c("STATx", "AP1"), verbose=TRUE)
}
\keyword{misc}
