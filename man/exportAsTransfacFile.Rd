
\name{exportAsTransfacFile}
\alias{exportAsTransfacFile}
\title{Write Transfac Files} 

\description{
Export an object of class motiv as a Transfac file type.
} 
\usage{
exportAsTransfacFile(x, file)
}
\arguments{
	\item{x}{An object of class \code{motiv} to be export.}
	\item{file}{A character string naming a file.}
}
\details{This function is made to provided standard output file used by STAMP. It take an object of class \code{motiv} and write two files named \code{*_matched.txt} and \code{*_match_pairs.txt} containing alignments and identified PWMs.

For more information about the Transfac file format, please refere to http://www.benoslab.pitt.edu/stamp/help.html.
}
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

#####Export#####
#exportAsTransfacFile(stat1.filter.combine, file="stat1_analysis")
}
\keyword{misc}
 
