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
	\S4method{plot}{motiv,ANY}(x, y=NULL, main=NULL, sub=NULL, ncol=0, nrow=0, top=3, bysim=TRUE, rev=FALSE,...)
	\S4method{plot}{motiv,gadem}(x, y, sort=FALSE, group=FALSE, main=NULL, sub=NULL, ncol=0, nrow=0, xlim=NULL, correction=TRUE, method=3, bysim=TRUE, strand=FALSE, bw="nrd0",  type="distribution",...)
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
	\item{bw}{The smoothing bandwidth to be used. (see ?density)}
	\item{xlim}{numeric vectors of length 2, giving the x coordinates ranges.}
	\item{method}{The method that should be used to calculate the motifs distances.}
	\item{correction}{If \code{TRUE}, corrects the position according to the alignment.}
	\item{...}{Further potential arguments passed to methods.}
}

\details{
	A single \code{motiv} object (usualy provied by \code{motifMatch}) will plot the list of identified transcription factors for each motif. 
With \code{rev=TRUE}, the transcription factor logo will be print to correspond to the real alignment instead of original TF PWM.
	
	Giving a \code{motiv} object and a \code{gadem} object with \code{type="distribution"} will show the motif repartition within gadem peaks.
If \code{strand=TRUE}, a distinct distribution is made for forward and reverse strand.

	A \code{var.test} is automatically made to help to distinguish centered distribution. The distribution with lowest variance is assign as "reference" distribution to compute the \code{var.test} statistic.
With \code{sort=TRUE}, distribution are plot according decreasing statistic.	

	\code{type="distance"} indicates to compute and plot the distance between each pair of motif. It aslo provied Venn diagramm that returns the proportion of common sequences per pair of motif. There are few methods available. It doesn't change for a wide however.
\code{method}=
\describe{
		\item{1}{Every motifs of a sequence are used to compute the distance between motif 1 and 2}
		\item{2}{Same as 1 but it return the distance of the motif 2 versus motif 1}
		\item{3}{Apply the method 4 or 5 depending of the lowest standard deviation}
		\item{4}{Get the position of the motif 1 in the middle and computes relative distance}
		\item{5}{Same as 4 but with the motif 2 centered.}
		}	
	The \code{group} argument indicates to consider similar motif as a single motif.
	
	With \code{correction=TRUE} the motif position is corrected accoring to the alignment. It means that the gap/"N" contained in the alignments are removed to give a corrected start and end position. 
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

#####Plots#####
plot(stat1.filter.combine, ncol=2,top=5, rev=FALSE, main="STAT1", bysim=TRUE)
plot(stat1.filter.combine ,gadem,ncol=2, type="distribution", correction=TRUE, group=FALSE, bysim=TRUE, strand=FALSE, sort=TRUE, main="STAT1")
plot(stat1.filter.combine ,gadem,type="distance", correction=TRUE, group=TRUE, bysim=TRUE, main="STAT1", strand=FALSE, method=3, xlim=c(-100,100), bw=8)
}
\keyword{methods}

