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
      similar <- names(object@input)[names==names(names.factor)[i]]
      cat("\t  -Motifs", paste( similar[1:(names.factor[i]-1)],"and"), similar[names.factor[i]], "are similar\n")
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

#####filter#####

setMethod(
"summary",
"filter",
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

#####filter#####
setMethod(
"summary",
"filters",
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
setMethod(
"summary",
"list",
function(object) {	
  if(all(lapply(object, class)=="motiv"))
  {
    for (i in 1:length(object))
    {
      print(names(object)[i])
      summary(object[[i]])	
    }
  }
  else if (all(lapply(object, class)=="filter") )
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
setGeneric("similar", function(x) standardGeneric("similar"))

setMethod("similar",
"motiv",
function(x){
  sim <- NULL
  for (i in 1:length(x@bestMatch))
  {
    sim <- c(sim, x@bestMatch[[i]]@similar)
  }
  return(unique(sim))
})

setMethod("similar",
"list",
function(x){
  sim <- NULL
  if (length(x)!=0)
  {
    for (i in 1:length(x))
    {
      sim=c(sim, x[[i]]@similarMotif)
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
function(x, i, j=ANY, bysim=TRUE, ..., exact=TRUE, drop=FALSE){
  selected <- NULL
  selectedMatch <- NULL
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
          selectedNames <- x@bestMatch[[k]]@similar
        } 
        else 
        {
          selectedNames <- x@bestMatch[[k]]@name
        }
        if (!exact)
        {
          if (grepl(i[kk], selectedNames))
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
    selectedMatch <- c(selectedMatch, new("matches", name=x@bestMatch[[l]]@name, aligns=x@bestMatch[[l]]@aligns, similar=x@bestMatch[[l]]@similar, valid=x@bestMatch[[l]]@valid)	)  	
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
  filt <- f@filters
  filterdMotiv <- filterMotiv(x, filt, exact, verbose)
  return(filterdMotiv)
})

#####combine#####
setGeneric("combine", function(x, f, name=NULL, exact=FALSE, verbose=TRUE) standardGeneric("combine"))

setMethod("combine",
signature(x="motiv", f="list"),
function(x, f, name=NULL, exact=FALSE, verbose=TRUE )
{
  motiv <- x
  combinedMotiv <- combineMotiv(motiv, f, name, verbose, exact)
  return(combinedMotiv)
})

setMethod("combine",
signature(x="motiv", f="filters"),
function(x, f, name=NULL, exact=FALSE, verbose=TRUE)
{
  f <- list(f)
  combine(x, f, name=name, verbose=verbose, exact=exact)
})

#####split#####

setMethod("split",
signature(x="motiv", f="list"),
function(x, f, exact=FALSE, drop = FALSE, verbose=TRUE, ...)
{
  splitedMotiv <- splitMotiv (x, f, drop, exact, verbose)
  return(splitedMotiv)
})

setMethod("split",
signature(x="motiv", f="filters"),
function(x, f, exact=FALSE, drop=FALSE, verbose=TRUE, ...)
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
function(x, y=NULL, main=NULL, sub=NULL, ncol=0, nrow=0, top=3, bysim=TRUE, rev=FALSE,...)
{
  if (top > length(x@bestMatch[[1]]@aligns) || top <1)
  {
    stop(paste("Top must be inferior to "), length(x@bestMatch[[1]]@aligns))
  }
  if (class(x)!="motiv")
  {
    stop("motiv must be an object of class motiv.")
  }
  if(ncol==0 && nrow==0)
  {
    ncol=ceiling(sqrt(length(x@bestMatch)))
    nrow=ncol
  } 
  else if (ncol==0 && nrow!=0)
  {
    ncol=ceiling(length(x@bestMatch)/nrow)
  }
  else if ((ncol!=0 && nrow==0))
  {
    nrow=ceiling(length(x@bestMatch)/ncol)
  }
  if ( (nrow*ncol) < length(x@bestMatch) )
  {
    print(paste("Warning : layout is not sufficient, some results couldn't be plot."))
  }
  plotMotiv(x, ncol, nrow, top, bysim, rev, main, sub)
})

#####motiv, gadem#####

setMethod(
"plot",
signature(x="motiv", y="gadem"),
function(x, y, sort=FALSE, group=FALSE, main=NULL, sub=NULL, ncol=0, nrow=0, xlim=NULL, correction=TRUE, method=3, bysim=TRUE, strand=FALSE, bw="nrd0",  type="distribution",...)
{
  if(type=="distance")
  {
    pos <- calculatePositionVector(x, y, group, correction)
    meanLength = y@motifList[[1]]@alignList[[1]]@end-y@motifList[[1]]@alignList[[1]]@start
    plotDistance( pos, strand, main, sub, FALSE, bysim, xlim, method, bw, meanLength)
  }
  else
  {
    for (g in 1:length(y@motifList))
    {
      motifs <- NULL
      if (y@motifList[[g]]@name %in% names(x@input))
      { motifs <- c(motifs, y@motifList[[g]])}
    }
    pos <- calculatePositionVector(x, y, group, correction)
    meanLength = y@motifList[[1]]@alignList[[1]]@end-y@motifList[[1]]@alignList[[1]]@start
    plotDistribution(pos, group, main, sort, ncol, nrow, strand, bysim, meanLength)
  }
})
