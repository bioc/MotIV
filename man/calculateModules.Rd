 
 
\name{calculateModules}
\alias{calculateModules}
\title{Calculate Modules} 

\description{
This function is use to calculate modules from two RangedData objects.} 
\usage{
  calculateModules(x, y, range=250, merge=TRUE)
}    
\arguments{
	\item{x}{A \code{RangedData} object.}
	\item{y}{A \code{RangedData} object.}
	\item{range}{The maximum distance allowed for a module.}
	\item{merge}{Logical. If \code{TRUE}, the two motifs composing the modules are merged.}
	}
\details{\code{findModule} calculates the close positions ie. modules between two RangedData.

If 'merge' is \code{TRUE}, the two motifs composing the modules are merged. It means gap between the two motifs are combined into the new module. If you only want to select true motifs, set 'merge=FALSE'.
} 
\author{Eloi Mercier <\email{eloi.mercier@ircm.qc.ca}>}
\value{A list of matrix.}
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
modules.stat1.ap1 <- calculateModules(stat1.rd, ap1.rd, range=250, merge=TRUE)
}
\keyword{misc}
 
