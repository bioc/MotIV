\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,motiv,ANY-method}
\alias{plot,motiv-method}
\alias{plot,motiv,gadem-method}
\alias{plot,motiv,gadem-method}
\alias{plot}


\title{Plot Motiv} 
\description{
This functions are used to vizualise and validate \code{motiv} analysis.
} 

\usage{
	\S4method{plot}{motiv,ANY}(x, y=NULL, main=NULL, sub=NULL, ncol=0, nrow=0, top=3, bysim=TRUE, rev=FALSE, trim=0.05, cex)
	\S4method{plot}{motiv,gadem}(x, y, sort=FALSE, group=FALSE, main=NULL, sub=NULL, ncol=0, nrow=0, xlim=NULL, correction=TRUE,  bysim=TRUE, strand=FALSE, type="distribution", trim=0.05, col=c("blue", "red"), border=c("black", "black"), lwd=2, lty=1, nclass=20, bw="nrd0", cex=1, vcol=c("red", "green"))
}

\arguments{
	\item{x}{An object of class \code{motiv}.}
	\item{y}{The GADEM type object associated with the \code{motiv} object.}
	\item{ncol, nrow}{A numeric value giving the the number of columns and rows to plot.}
	\item{top}{A numeric value giving the number of best matches per motif to display.}
	\item{rev}{A logical value. If \code{TRUE}, print reverse motif for negatif strand.}	
	\item{main}{An overall title for the plot: see \code{title}.}
	\item{sub}{A sub title for the plot: see 'title'}	
	\item{type}{What type of plot should be drawn. Possible values are : distribution to display the binding sites distribution within the peaks or distance to show the pairwise distance between motifs.}
	\item{strand}{If \code{TRUE}, distribution will be plot for both forward and reverse strand.}
	\item{group}{If \code{TRUE}, similar motifs will be grouped.}
	\item{sort}{If \code{TRUE}, motifs will be plot according their computed variance.}
	\item{bysim}{If \code{TRUE}, the 'similar' field (defined with the \code{combine} function) will be print instead of the original name.}
	\item{xlim}{A numeric vectors of length 2, giving the x coordinates ranges.}
	\item{correction}{If \code{TRUE}, corrects the position according to the alignment.}
	\item{trim}{A \code{numeric} value. Define the mimimun information content value for which the logo letters are shown.}
	\item{col, border, lwd, lty}{Define respectively the color, the border, the line wide and the line type of both curve and histogram. See 'par'.}
	\item{nclass}{A numerical value giving the number of class for the histogram.}
	\item{bw}{he smoothing bandwidth to be used to calculate the density. See \code{density}.}
	\item{cex, vcol}{A numerical value giving the amount by which plotting text should be magnified relative to the default.}
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

#####Plots#####
plot(foxa1.filter.combine, ncol=2,top=5, rev=FALSE, main="FOXA", bysim=TRUE)
plot(foxa1.filter.combine ,gadem,ncol=2, type="distribution", correction=TRUE, group=FALSE, bysim=TRUE, strand=FALSE, sort=TRUE, main="FOXA", nclass=20, bw=2)
plot(foxa1.filter.combine ,gadem,type="distance", correction=TRUE, group=TRUE, bysim=TRUE, main="FOXA", strand=FALSE, xlim=c(-100,100), nclass=20, bw=8) 
}
\keyword{methods}

