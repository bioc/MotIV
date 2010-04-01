getSeqID <- function(x) 
{
	seqID <- NULL
		for (i in 1:length(x@alignList))
			{seqID <- c(seqID , x@alignList[[i]]@seqID)}
	return(seqID)
}

occurences <- function(gadem) 
{
	modulesID <- lapply(gadem@motifList,getSeqID)	
	table.modules <- list()
	for (i in seq(length(modulesID)))
	{
		table.modules[[i]]<-sapply( seq(max(unlist(modulesID))),function(x) length(modulesID[[i]][modulesID[[i]]==x]))
	}
	table <- data.frame(table.modules)
	names(table) <- names(gadem)
	return(table)
}

cooccurences <- function(x)
{
	table=matrix(ncol=dim(x)[2], nrow=dim(x)[2])
	colnames(table) <- names(x)
	rownames(table) <- names(x)
	for (j in seq(dim(x)[2]))
	{
		for (i in seq(dim(x)[2]))
		{
			table[j,i]=dim(x[x[,i]!=0&x[,j]!=0,])[1]
		}
	}
	return(table)
}