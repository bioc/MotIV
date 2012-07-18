\name{seqLogo2}
\alias{seqLogo2}
\title{Plot a sequence logo for a given position weight matrix}
\description{
  This function takes the 4xW position weight matrix of a DNA sequence
  motif and plots the corresponding sequence logo. 
}
\usage{
  seqLogo2(pwm, ic.scale=TRUE, xaxis=TRUE, yaxis=TRUE, xfontsize=15, yfontsize=15,vmargins=c(0,0),hmargins=c(0,0),size=1, trim=0)
}
\arguments{
	\item{pwm}{\code{numeric} The 4xW position weight matrix.}
	\item{ic.scale}{\code{logical} If TRUE, the height of each column is proportional to its information content. Otherwise, all columns have the same height.}
	\item{xaxis}{\code{logical} If TRUE, an X-axis will be plotted.}
	\item{yaxis}{\code{logical} If TRUE, a Y-axis will be plotted.}
	\item{xfontsize}{\code{numeric} Font size to be used for the X-axis.}
	\item{yfontsize}{\code{numeric} Font size to be used for the Y-axis.} 
	\item{vmargins}{\code{numeric} Vertical margins.} 
	\item{hmargins}{\code{numeric} Horizontal margins.} 
	\item{size}{\code{numeric} Graphic size.} 
	\item{trim}{\code{numeric} Print nucleotide only if the infomation content is superior to this \code{trim} threshold.}
}
\details{
  Within each column, the height of a given letter is proportional to its frequency at that position.
  If ic.scale is TRUE, the height of each column in the plot indicates the information content at that position of the motif. Otherwise, the height of all columns are identical.
  
  This is an internal function for the package MotIV. User should prefer the \code{seqLogo} function from the package \code{seqLogo} to vizualise individual motif.
}
\value{
  None.
}
\author{Oliver Bembom, \email{bembom@berkeley.edu}, modified by Eloi Mercier <\email{emercier@chibi.ubc.ca}>}

\examples{
#mFile <- system.file("Exfiles/pwm1", package="seqLogo")
#m <- read.table(mFile)
#pwm <- makePWM(m)
#seqLogo2(pwm)
}
\keyword{misc}
