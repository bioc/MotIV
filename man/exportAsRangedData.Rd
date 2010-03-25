
\name{exportAsRangedData}
\alias{exportAsRangedData}
\title{Export MotIV Results} 

\description{
Export your 
} 
\usage{
exportAsRangedData(x, y, bysim=TRUE, correction=TRUE)
}
\arguments{
	\item{x}{An object of class \code{motiv}.}
	\item{y}{The rGADEM type object associated with the \code{motiv }object.}
	\item{bysim}{If \code{TRUE}, the 'similar' field (defined with the \code{combine} function) will be print instead of the original name.}
	\item{correction}{If \code{TRUE}, corrects the position according to the alignment.}	
}
\value{An object of type \code{RangedData}.}
\details{Use this function to export the results into a \code{RangedData} object.}
\author{Eloi Mercier <\email{eloi.mercier@ircm.qc.ca}>}
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

#####Plots#####
#plot(stat1.filter.combine, ncol=2,top=5, rev=FALSE, main="STAT1", bysim=TRUE)
#plot(stat1.filter.combine ,gadem,ncol=2, type="distribution", correction=TRUE, group=FALSE, bysim=TRUE, strand=FALSE, sort=TRUE, main="STAT1")
#plot(stat1.filter.combine ,gadem,type="distance", correction=TRUE, group=TRUE, bysim=TRUE, main="STAT1", strand=FALSE, method=3, xlim=c(-100,100), bw=8)

#####RangedData#####
stat1.rd <- exportAsRangedData(stat1.filter.combine["STATx"], gadem)
ap1.rd <- exportAsRangedData(stat1.filter.combine["AP1"], gadem)
#modules.stat1.ap1 <- calculateModules(stat1.rd, ap1.rd, range=250, merge=TRUE)
}
\keyword{misc}
 
