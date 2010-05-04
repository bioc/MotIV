
#####RANGED DATA#####
exportAsRangedData <- function (x, y, bysim=TRUE, correction=TRUE)
{
	if (class(x)!="motiv")
	{
		stop("x must be an object of class motiv.")
	}
	vector.pos <- calculatePositionVector(x, y, group=F, correction=correction)
	curr <- which(names(y)%in%names(x))
	
	data <- RangedData()
	
	p=1
	for (i in curr)
	{
		data <- rbind(data, vector.pos[[p]]@positionVector[,"strand"])	
		p=p+1
	}	
	return(data)
}
