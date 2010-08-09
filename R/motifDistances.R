
motifDistances <- function (inputPWM, DBscores=jaspar.scores, cc="PCC", align="SWU", top=5, go=1, ge=0.5) {
	res<-.Call("RmotifDistances", cc=cc, align=align, top=top, go=go, ge=ge, inputPWM=inputPWM, inputScores=DBscores)
	if (!is.null(res))
	{
		colnames(res)=names(inputPWM)
		rownames(res)=names(inputPWM)
		return(as.dist(res))
	}
}

###HCLUST###

motifHclust <- function(x,...){
hc <- do.call("hclust", list(x,...))
	return(hc)
}

####CUTREE###

motifCutree <- function(tree,k=NULL, h=NULL){
	cut<-do.call("cutree", list(tree=tree, k=k, h=h))
	f.cut <- list()
	for (i in unique(cut))
	{
		f.cut[[i]]<- setFilter(name=names(cut[cut==i]))
	}
	return(f.cut)
}

