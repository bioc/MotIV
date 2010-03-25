\name{transcriptionFactor-class}
\docType{class}
\alias{transcriptionFactor-class}

\title{Transcription Factor Class}
\description{This object contains the Transcription Factor name and PWM.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("transcriptionFactor", name, pwm)}.
}


\section{Slots}{
\describe{
	\item{name}{TF name.}
	\item{pwm}{TF PWM.}
}}


\author{Eloi Mercier <\email{eloi.mercier@ircm.qc.ca}>}
\seealso{
	 \code{motiv} , \code{matches}, \code{alignments}
}
\examples{
showClass("transcriptionFactor")
}
\keyword{classes}
