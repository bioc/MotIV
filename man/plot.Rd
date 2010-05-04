\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,motiv,ANY-method}
\alias{plot,motiv-method}
\alias{plot,motiv,gadem-method}
\alias{plot,motiv,gadem-method}

\title{Plot Motiv} 
\description{
This functions are used to vizualise and validate \code{motiv} analysis.
} 

\usage{
	\S4method{plot}{motiv,ANY}(x, y=NULL, main=NULL, sub=NULL, ncol=0, nrow=0, top=3, bysim=TRUE, rev=FALSE, trim=0.05, ...)
	\S4method{plot}{motiv,gadem}(x, y, sort=FALSE, group=FALSE, main=NULL, sub=NULL, ncol=0, nrow=0, xlim=NULL, correction=TRUE,  bysim=TRUE, strand=FALSE, type="distribution", trim=0.05, harg=list(), darg=list(), carg=list(), varg=list(), ...)
}

\arguments{
	\item{x}{An object of class \code{motiv}.}
	\item{y}{The GADEM type object associated with the \code{motiv} object.}
	\item{ncol}{The desired number of columns.}
	\item{nrow}{The desired number of rows.}	
	\item{top}{The top matches to be print.}
	\item{rev}{Print reverse motif for negatif strand.}	
	\item{main}{A title of the plot.}
	\item{sub}{A sub title of the plot.}	
	\item{type}{The type of plotting.}
	\item{strand}{If TRUE, two distributions will be plot for forward and reverse strand.}
	\item{group}{If TRUE, similar motif will be group.}
	\item{sort}{If TRUE, motifs will be plot by decreased variance order.}
	\item{bysim}{If \code{TRUE}, the 'similar' field (defined with the \code{combine} function) will be print instead of the original name.}
	\item{xlim}{numeric vectors of length 2, giving the x coordinates ranges.}
	\item{correction}{If \code{TRUE}, corrects the position according to the alignment.}
	\item{harg, darg, carg, varg}{Parameters for histograms, densities, curves and venn diagrams.}
	\item{trim}{\code{numeric} Only print letters with information content superior to a treshold given by \code{trim}}
	\item{...}{Further potential arguments passed to methods.}
}

\details{
	A single \code{motiv} object (usualy provied by \code{motifMatch}) will plot the list of identified transcription factors for each motif. 
With \code{rev=TRUE}, the transcription factor logo will be print to correspond to the real alignment instead of original TF PWM.
	
	Giving a \code{motiv} object and a \code{gadem} object with \code{type="distribution"} will show the motif repartition within gadem peaks.
If \code{strand=TRUE}, a distinct distribution is made for forward and reverse strand.

	A \code{var.test} is automatically made to help to distinguish centered distribution. The distribution with lowest variance is assign as "reference" distribution to compute the \code{var.test} statistic.
With \code{sort=TRUE}, distribution are plot according decreasing statistic.	

	\code{type="distance"} indicates to compute and plot the distance between each pair of motif. It aslo provied Venn diagramm that returns the proportion of common sequences per pair of motif.

	The \code{group} argument indicates to consider similar motif as a single motif.
	
	With \code{correction=TRUE} the motif position is corrected accoring to the alignment. It means that the gap/"N" contained in the alignments are removed to give a corrected start and end position. 
}
\author{Eloi Mercier <\email{eloi.mercier@ircm.qc.ca}>}

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
foxa1.filter.combine <- combine(foxa1.filter, c(f.foxa1, f.ap1), exact=FALSE, name=c("FOXA1", "AP1"), verbose=TRUE)

#####Plots#####
plot(foxa1.filter.combine, ncol=2,top=5, rev=FALSE, main="FOXA", bysim=TRUE)
plot(foxa1.filter.combine ,gadem,ncol=2, type="distribution", correction=TRUE, group=FALSE, bysim=TRUE, strand=FALSE, sort=TRUE, main="FOXA", harg=list(nint=20), carg=list(bw=2))
plot(foxa1.filter.combine ,gadem,type="distance", correction=TRUE, group=TRUE, bysim=TRUE, main="FOXA", strand=FALSE, xlim=c(-100,100), harg=list(nint=20), darg=list(bw=8), carg=list(bw=2), varg=list(col=c("red", "blue")))
}
\keyword{methods}

