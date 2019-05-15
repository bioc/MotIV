
#####GRanges#####
exportAsGRanges <- function (x, y, correction=TRUE)
{
	if (class(x)!="motiv")
	{
		stop("x must be an object of class motiv.")
	}
	vector.pos <- calculatePositionVector(x, y, group=F, correction=correction)
	curr <- which(names(y)%in%names(x))
	
	data <- GRanges()
	
	p=1
	for (i in curr)
	{
		data <- c(data, granges(vector.pos[[p]]@positionVector))	
		p=p+1
	}	
	return(data)
}
