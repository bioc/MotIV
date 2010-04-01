
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
		data.ir=IRanges(start=sapply(y@motifList[[i]]@alignList, function(x){x@start})+vector.pos[[p]]@positionVector$start, end=sapply(y@motifList[[i]]@alignList, function(x){x@start})+vector.pos[[p]]@positionVector$end)
		chr <- sapply(y@motifList[[i]]@alignList, function(x){x@chr})
		data = rbind(data,RangedData(data.ir, space=gsub("chr", "", chr), strand=sapply(y@motifList[[i]]@alignList, function(x){x@strand})))
		p=p+1
	}	
	return(data)
}
