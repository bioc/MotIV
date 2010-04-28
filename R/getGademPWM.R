
getGademPWM <- function (y)
{
 .Deprecated("getPWM", package="rGADEM")
	PWM <- sapply(y@motifList, 
	function(x){
		pwm=x@pwm
		colnames(pwm) <- 1:(length(pwm)/4)
		rownames(pwm) <- c("A","C","G","T")
		pwm
	})  
	names(PWM)<- names(y)
	return(PWM)
}




