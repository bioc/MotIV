
\name{exportAsTransfacFile}
\alias{exportAsTransfacFile}
\alias{exportAsTransfacFile-methods}
\alias{exportAsTransfacFile,motiv-method}
\alias{exportAsTransfacFile,list-method}
\title{Write Transfac Files} 

\description{
Export an object of class motiv as a Transfac file type.
} 
\usage{
	\S4method{exportAsTransfacFile}{motiv}(x, file)
	\S4method{exportAsTransfacFile}{list}(x, file)
}
\arguments{
	\item{x}{An object of class \code{motiv} to be export.}
	\item{file}{A character string naming a file.}
}
\details{
This function is made to provide standard output file used by STAMP. It take an object of class \code{motiv} and write two files named \code{*_matched.txt} and \code{*_match_pairs.txt} containing alignments and identified PWMs.

For more information about the Transfac file format, please refere to http://www.benoslab.pitt.edu/stamp/help.html.
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
foxa1.filter.combine <- combineMotifs(foxa1.filter, c(f.foxa1, f.ap1), exact=FALSE, name=c("FOXA1", "AP1"), verbose=TRUE)

#####Export#####
#exportAsTransfacFile(foxa1.filter.combine, file="foxa1_analysis")
}
\keyword{misc}
 
