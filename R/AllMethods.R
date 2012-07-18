###############################
##########SUMMARY###############
###############################

#####motiv#####
setMethod("summary", signature("motiv"),
function(object) {
	argv <- object@argv
	input.names <- names(object@input)
	motifs <- length(object@bestMatch)
	bestMatch <- length(object@bestMatch[[1]]@aligns)
	factor.names <- NULL
	names <- NULL
	for (i in 1:motifs)
	{
		for (j in 1:bestMatch)
		{
			factor.names <- c(factor.names, object@bestMatch[[i]]@aligns[[j]]@TF@name)
		}
		names <- c(names, object@bestMatch[[i]]@name)
	}
	cat("\tNumber of input motifs : ", motifs, "\n")
	cat("\tInput motifs names : ", input.names, "\n")
	names.factor <- summary(factor(names))
	for (i in 1:length(names.factor))
	{
		if (names.factor[i]>1)
		{
			similarity <- names(object@input)[names==names(names.factor)[i]]
			cat("\t  -Motifs", paste( similarity[1:(names.factor[i]-1)],"and"), similarity[names.factor[i]], "are similar\n")
		}
	}
	cat("\tNumber of matches per motif: ", bestMatch, "\n")	
	cat("\tMatches repartition : \n")
	print(sort(summary(factor(factor.names)), decreasing=T))
	cat("\tArguments used : \n")
	cat("\t  -metric name : ", argv["cc"], "\n")
	cat("\t  -alignment : ", argv["align"], "\n")
	if (argv["align"] != "SWU")
	{
		cat("\t  -gap open penality : ", argv["go"], "\n")
		cat("\t  -gap extension penality : ", argv["ge"], "\n")
	}
	cat("\n")
})

#as.data.frame
setMethod("as.data.frame", signature("motiv"),
function (x) {
	if (length (x@bestMatch) == 0)
	return (NA)
	df = data.frame (stringsAsFactors=FALSE)
	for (k in seq(length(x@bestMatch)))
	{
		alignments = x@bestMatch[[k]]@aligns
		motif = x@bestMatch[[k]]@name
		for (alignment in alignments) {
			TF = alignment@TF@name
			eVal = alignment@evalue
			sequence = alignment@sequence
			match = alignment@match
			strand = alignment@strand
			df = rbind (df, data.frame ( motif=motif, TF=match, eVal=eVal, sequence=sequence, match=match, strand=strand, stringsAsFactors=FALSE))
		}        
	}
	return (df)
})

#####filter#####
setMethod("summary", "filter",
function(object){
	for (i in 1:length(object@name))
	{
		if(length(unlist(object@name[[i]]))>1)
		{
			cat("\tNames = ", object@name[[i]][1], paste( "or", object@name[[i]][-1]),"\n")
		} 
		else if(object@name[[i]]!="")
		{
			cat("\tName = ", unlist(object@name),"\n")
		}
		if(object@tfname[[i]]!="")
		{cat("\tTF name = ", object@tfname[[i]],"\n")}
		if(object@evalueMax[[i]] < 1)
		{cat("\tEvalue max = ", object@evalueMax[[i]],"\n")}
		if(object@lengthMax[[i]] != 100)
		{cat("\tLength max = ", object@lengthMax[[i]],"\n")}
		if (unlist(object@top[[i]]) < 10)
		{cat("\tin the ", object@top[[i]], " first positions.\n")}
		if (i<length(object@name))
		{cat("\t\tAND\n")}
	}
})

#####filters#####
setMethod("summary", "filters",
function(object){
	cat("Filter will select motifs that satisfy the conditions :\n")
	for(i in 1:length(object@filters))
	{
		summary(object@filters[[i]])
		if(i<length(object@filters))
		{
			cat("OR\n")
		}
	}
})

#####list#####
setMethod("summary", "list",
function(object) {	
	if(all(lapply(object, class)=="motiv"))
	{
		for (i in 1:length(object))
		{
			print(names(object)[i])
			summary(object[[i]])	
		}
	}
	else if (all(lapply(object, class)=="filters") )
	{
		for (i in 1:length(object))
		{
			summary(object[[i]])
		}
	}
})

#####motifs#####
setGeneric("viewMotifs", function(x, n=100) standardGeneric("viewMotifs"))
setMethod("viewMotifs",
"motiv",
function(x, n=100){
	factor.names <- NULL
	motifs <- length(x@bestMatch)
	bestMatch <- length(x@bestMatch[[1]]@aligns)
	for (i in 1:motifs)
	{
		for (j in 1:bestMatch)
		{
			factor.names <- c(factor.names, x@bestMatch[[i]]@aligns[[j]]@TF@name)}
    }
    return(sort(summary(factor(factor.names), n+1), decreasing=T))
})

###############################
###########NAMES###############
###############################

#####acces#####
setMethod("names",
"motiv",
function(x){
	names<-sapply(x@bestMatch,function(x){x@name})
	return(names)
})

setMethod("names",
"filter",
function(x){
	return(x@name[[1]])
})

setMethod("names",
"filters",
function(x){
	return(names(x@filters[[1]]))
})

#####replace#####
setReplaceMethod("names",
"motiv",
function(x, value){
	for (i in 1:length(x@bestMatch))
	{
		x@bestMatch[[i]]@name <- value[i]
	}
	return(x)
})

#####similarity#####
setGeneric("similarity", function(x) standardGeneric("similarity"))

setMethod("similarity",
"motiv",
function(x){
	sim <- NULL
	for (i in 1:length(x@bestMatch))
	{
		sim <- c(sim, x@bestMatch[[i]]@similarity)
	}
	return(sim)
})

setMethod("similarity",
"list",
function(x){
	sim <- NULL
	if (length(x)!=0)
	{
		for (i in 1:length(x))
		{
			sim=c(sim, x[[i]]@similarity)
		}	
	}
	return(sim)
})

###############################
###########LENGTH###############
###############################

setMethod("length",
"motiv",
function(x){
	res <- length(x@bestMatch)
	return(res)
})

###############################
###########SELECT###############
###############################

setMethod("[",
"motiv",
function(x, i, j=ANY, bysim=TRUE, ..., exact=TRUE, ignore.case=FALSE, drop=FALSE){
	selected <- NULL
	selectedMatch <- NULL
	if (ignore.case) exact <- FALSE
	if (is.numeric(i))
	{
		selected = i
	}
	else 
	{
		for (kk in 1:length(i))
		{
			for (k in 1:length(x@bestMatch))
			{
				if (bysim)
				{
					selectedNames <- x@bestMatch[[k]]@similarity
				} 
				else 
				{
					selectedNames <- x@bestMatch[[k]]@name
				}
				if (!exact)
				{
					if (grepl(i[kk], selectedNames,  ignore.case= ignore.case))
					{
						selected <- c(selected, k)
					}
				} 
				else 
				{
					if (i[kk] == selectedNames)
					{
						selected <- c(selected, k)
					}
				}
			}
		}
	}
	if (drop)
	{
		selected <- (1:length(x))[-selected]
	}
	for (l in unique(selected))
	{
		selectedMatch <- c(selectedMatch, new("matches", name=x@bestMatch[[l]]@name, aligns=x@bestMatch[[l]]@aligns, similarity=x@bestMatch[[l]]@similarity, valid=x@bestMatch[[l]]@valid)	)  	
	}		
	if (!is.null(selectedMatch))
	{
		res <- new("motiv", input=x@input[selected], bestMatch=selectedMatch, argv=x@argv)
	} 
	else 
	{
		res <- NULL
	}
	return(res)
})

###############################
############SHOW###############
###############################

#####motiv#####
setMethod("show", "motiv",
function(object)
{
	cat("\tObject of class 'motiv'","\n")
	cat("\tThis object has the following slots: \n")
	cat("\targv, bestMatch, input\n\n")
})

#####filter#####
setMethod("show", "filter",
function(object)
{
	cat("\tObject of class 'filter'","\n")
	cat("\tThis object has the following slots: \n")
	cat("\tname, tfname, top, evalueMax, lengthMax\n\n")
})

###############################
###########FILTERS###############
###############################

#####SetFilter#####
setFilter <- function (name="", tfname="", evalueMax=1, top=10, lengthMax=100, valid=NULL)
{
	if (name=="" && tfname=="" && evalueMax>=1 && top>=10 && lengthMax>=100)
	{
		warning("Filter does not contain condition.")
	}
	fit <- new("filter", name=list(name), tfname=list(tfname), evalueMax=list(evalueMax), top=list(top), lengthMax=list(lengthMax), valid=list(valid))
	filters <- new("filters", filters=list(fit))
	return(filters)
}

#####&#####
setMethod("&",
signature(e1="filters", e2="filters"),
function(e1, e2){
	fit <- list()
	p=1
	for(i in 1:length(e1@filters))
	{
		for (j in 1:length(e2@filters))
		{
			fit[[p]] <- new("filter", name=c(e1@filters[[i]]@name, e2@filters[[j]]@name), tfname=c(e1@filters[[i]]@tfname, e2@filters[[j]]@tfname), top=c(e1@filters[[i]]@top, e2@filters[[j]]@top), evalueMax=c(e1@filters[[i]]@evalueMax, e2@filters[[j]]@evalueMax), lengthMax=c(e1@filters[[i]]@lengthMax, e2@filters[[j]]@lengthMax),
							valid=c(e1@filters[[i]]@valid, e2@filters[[j]]@valid))
			p=p+1
		}
	}
	filters <- new("filters", filters=fit)
	return(filters)
})

#####|#####

setMethod("|",
signature(e1="filters", e2="filters"),
function(e1, e2){
	fit <- c(e1@filters, e2@filters)
	filters <- new("filters", filters=fit)
	return(filters)
})

#####filter#####

setGeneric("filter", function(x, f, exact=FALSE, verbose=TRUE) standardGeneric("filter"))

setMethod(
"filter",
signature(x="motiv", f="filters"),
function(x, f, exact=FALSE, verbose=TRUE)
{
	filteredMotiv <- filterMotiv(x, f@filters, exact, verbose)
	return(filteredMotiv)
})

setMethod(
"filter",
signature(x="motiv", f="list"),
function(x, f, exact=FALSE, verbose=TRUE)
{
	ff=f[[1]]
	if (length(f) >=2)
{
	for (i in 2:length(f)){
	ff=ff|f[[i]]}
}
	filteredMotiv <- filterMotiv(x, ff@filters, exact, verbose)
	return(filteredMotiv)
})

#####combine#####
setGeneric("combineMotifs", function(x, y, name=NULL, exact=TRUE, verbose=TRUE) standardGeneric("combineMotifs"))

setMethod("combineMotifs",
signature(x="motiv", y="list"),
function(x, y, name=NULL, exact=TRUE, verbose=TRUE)
{
	motiv <- x
	combinedMotiv <- combineMotiv(motiv, y, name, verbose, exact)
	return(combinedMotiv)
})

setMethod("combineMotifs",
signature(x="motiv", y="filters"),
function(x, y, name=NULL, exact=TRUE, verbose=TRUE)
{
	y <- list(y)
	combineMotifs(x, y, name=name, verbose=verbose, exact=exact)
})

#####split#####

setMethod("split",
signature(x="motiv", f="list"),
function(x, f, exact=TRUE, drop = FALSE, verbose=TRUE, ...)
{
	splitedMotiv <- splitMotiv (x, f, drop, exact, verbose)
	return(splitedMotiv)
})

setMethod("split",
signature(x="motiv", f="filters"),
function(x, f, exact=TRUE, drop=FALSE, verbose=TRUE, ...)
{
	f <- list(f)
	split(x, f, drop=drop, exact=exact, verbose=verbose)
})

###############################
############PLOT################
###############################

#####motiv#####
setMethod(
"plot",
signature(x="motiv", y="ANY"),
function(x, y=NULL, main=NULL, sub=NULL, ncol=0, nrow=0, top=3, bysim=TRUE, rev=FALSE, trim=0.05, cex=1)
{
	plotMotiv(x, ncol, nrow, top, bysim, rev, main, sub, trim, cex)
})

#####motiv, gadem#####

setMethod(
"plot",
signature(x="motiv", y="gadem"),
function(x, y, sort=FALSE, group=FALSE, main=NULL, sub=NULL, ncol=0, nrow=0, xlim=NULL, correction=TRUE, bysim=TRUE, strand=FALSE,  type="distribution", trim=0.05, col=c("blue", "red"), border=c("black", "black"), lwd=2, lty=1, nclass=20, bw="nrd0", cex=1, vcol=c("red", "green"))
{
	pos <- calculatePositionVector(x, y,  group, correction)

	if(type=="distance")
	{
		table <- occurences(y[names(x)])[names(x)]
		
		if (group)
		{
			names(table) <- similarity(x)
			table <- as.data.frame(sapply(similarity(pos),function(x) {
				similar <- which(names(table)==x)
				if (length(similar)>1) {
				apply(table[,similar],1,sum)
				} else {
				table[,similar]
			}}))
		}
		plotDistance( pos, table, strand, main, FALSE, bysim, xlim, y@parameters[[1]]@nSequences, col, border, lwd, lty, nclass, bw, cex, vcol)
	}
	else
	{
		for (g in seq(nMotifs(y)))
		{
			motifs <- NULL
			if (y@motifList[[g]]@name %in% names(x@input))
			{ motifs <- c(motifs, y@motifList[[g]])}
		}
		plotDistribution(pos, group, main, sort, ncol, nrow, strand, bysim, trim, col, border, lwd, lty, nclass, bw, cex)
	}
})


###############################
###########EXPORT###############
###############################

setGeneric("exportAsTransfacFile", function(x, file) standardGeneric("exportAsTransfacFile"))

###PWMs###

setMethod("exportAsTransfacFile",
signature(x="list"),
function (x, file){ 
  transfac <- file (file, "w") 
  for (i in 1:length(x))
  {
      writeChar(paste("DE\t", names(x)[i],"\n", sep=""), transfac, eos=NULL)
      pwm <- format( round(x[[i]],4), digits=4)
      for (l in 1:dim(pwm)[2])
      {
        writeChar(paste(c(l-1, pwm[(4*l-4+1):(4*l)]), c("\t","\t","\t","\t","\n"), sep=""), transfac, eos=NULL)
      }
      writeChar("XX\n", transfac, eos=NULL)
  
  }
  close(transfac)
})

###motiv###

setMethod("exportAsTransfacFile",
signature(x="motiv"),
 function (x, file)
{ 
  if (class(x)!="motiv")
  {stop("motiv must be an object of class motiv.")}

  transfac <- file (paste(file, "_matched.txt", sep=""), "w") #PWMs found
  pairs <-   file (paste(file, "_match_pairs.txt", sep=""), "w") #alignments

  for (i in 1:length(x@input))
  {
    writeChar(paste(">\t", names(x@input)[i], "\n", sep=""), pairs, eos=NULL)
    for (j in 1:length(x@bestMatch[[i]]@aligns))
    {
      writeChar(paste(x@bestMatch[[i]]@aligns[[j]]@TF@name, "\t", format(x@bestMatch[[i]]@aligns[[j]]@evalue, digits=5, scientific=T), "\t", x@bestMatch[[i]]@aligns[[j]]@sequence, "\t", x@bestMatch[[i]]@aligns[[j]]@match, "\n", sep=""), pairs, eos=NULL)
      writeChar(paste("DE\t", x@bestMatch[[i]]@aligns[[j]]@TF@name,"\n", sep=""), transfac, eos=NULL)
      pwm <- format( round(x@bestMatch[[i]]@aligns[[j]]@TF@pwm,4), digits=4)
      for (l in 1:dim(x@bestMatch[[i]]@aligns[[j]]@TF@pwm)[2])
      {
        writeChar(paste(c(l-1, pwm[(4*l-4+1):(4*l)]), c("\t","\t","\t","\t","\n"), sep=""), transfac, eos=NULL)
      }
      writeChar("XX\n", transfac, eos=NULL)
    }
  }
  close(transfac)
  close(pairs)
  print(paste(file,"_matched_transfac.txt created.", sep=""))
  print(paste(file,"_match_pairs.txt created.", sep=""))
})

###############################
###########GETPWM###############
###############################

setMethod("getPWM",
signature(x="motiv"),
function (x){ 
  return(x@input)
})
