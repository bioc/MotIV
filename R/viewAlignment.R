

viewAlignments <- function (x)
{
	align <- list()
	for (i in 1:length(x@bestMatch))
	{
		seq <-{}
		match <-{}
		name <-{}
		evalue <-{}
		for (j in 1:length(x@bestMatch[[i]]@aligns))
		{
			seq <- c(seq, x@bestMatch[[i]]@aligns[[j]]@sequence)
			match <- c(match, x@bestMatch[[i]]@aligns[[j]]@match)
			name <- c(name, x@bestMatch[[i]]@aligns[[j]]@TF@name)
			evalue <- c(evalue, format(x@bestMatch[[i]]@aligns[[j]]@evalue, scientific=T, digits=5))
		}
	align[[i]]<- rbind(seq, match, evalue)
	names(align)[i]<- x@bestMatch[[i]]@name
	colnames(align[[i]])<- name
	}
	return(noquote(align))
} 
