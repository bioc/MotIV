
\name{exportAsGRanges}
\alias{exportAsGRanges}
\title{Export MotIV Results} 

\description{
Export your 
} 
\usage{
exportAsGRanges(x, y, correction=TRUE)
}
\arguments{
	\item{x}{An object of class \code{motiv}.}
	\item{y}{The rGADEM type object associated with the \code{motiv }object.}
	\item{correction}{If \code{TRUE}, corrects the position according to the alignment.}	
}
\value{An object of type \code{GRanges}.}
\details{Use this function to export the results into a \code{GRanges} object.}
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
summary(foxa1.analysis.jaspar )

#####Filters#####
f.foxa1<-setFilter(name="", tfname="FOXA1", top=3, evalueMax=10^-5)
f.ap1 <- setFilter (tfname="AP1", top=3)
f.foxa1.ap1 <- f.foxa1 | f.ap1 
foxa1.filter <- filter(foxa1.analysis.jaspar, f.foxa1.ap1, exact=FALSE, verbose=TRUE)
foxa1.split <- split(foxa1.analysis.jaspar, c(f.foxa1, f.ap1) , drop=FALSE, exact=FALSE, verbose=TRUE)
foxa1.filter.combine <- combineMotifs(foxa1.filter, c(f.foxa1, f.ap1), exact=FALSE, name=c("FOXA1", "AP1"), verbose=TRUE)

#####Plots#####
#plot(foxa1.filter.combine, ncol=2,top=5, rev=FALSE, main="FOXA1", bysim=TRUE)
#plot(foxa1.filter.combine ,gadem,ncol=2, type="distribution", correction=TRUE, group=FALSE, bysim=TRUE, strand=FALSE, sort=TRUE, main="FOXA1")
#plot(foxa1.filter.combine ,gadem,type="distance", correction=TRUE, group=TRUE, bysim=TRUE, main="FOXA1", strand=FALSE, method=3, xlim=c(-100,100), darg=list(bw=8))

#####GRanges#####
foxa1.gr <- exportAsGRanges(foxa1.filter.combine["FOXA1"], gadem)
ap1.gr <- exportAsGRanges(foxa1.filter.combine["AP1"], gadem)
}
\keyword{misc}
 
