
#####FILTER#####

filterMotiv <- function (x, f, exact, verbose)
{
	length.motif=length(x@bestMatch)
	selectedMotifs <- NULL
	selectedMatch <- NULL
	countSelectedMotifs=1	
	validMotifs <- NULL
	
	if (verbose)
	{
		message("motiv object contains ", length.motif, " motifs.")
	}
	
	for (filterNumber in 1:length(f))
	{
		countvalidMotifs=1
		select <- list()
		for (h in 1:length(f[[filterNumber]]@name))
		{
			select[[h]]=-1			
			for (i in 1:length.motif)
			{
				name <- FALSE
				TFname <- FALSE
				eval <- FALSE
				topX <- FALSE
				len <- TRUE
				
#check motif name
				if (exact && f[[filterNumber]]@name[[h]]!="")
				{
					isInNames <- any(unlist(lapply(f[[filterNumber]]@name[[h]], function(z){z==x@bestMatch[[i]]@name})))
				} 
				else 
				{
					isInNames <- any(unlist(lapply(f[[filterNumber]]@name[[h]], function(z){grepl(z, x@bestMatch[[i]]@name)})))
				}
				if (isInNames)
				{
					name <- TRUE	
				}
				
				validMotifs.tmp=NULL
				for (j in 1:length(x@bestMatch[[i]]@aligns))
				{			
#check TF names
					if (exact && f[[filterNumber]]@tfname[[h]]!="")
					{
						isInTFNames <- any(unlist(lapply(f[[filterNumber]]@tfname[[h]], function(z){z==x@bestMatch[[i]]@aligns[[j]]@TF@name})))
					} 
					else
					{
						isInTFNames <- any(unlist(lapply(f[[filterNumber]]@tfname[[h]], function(z){grepl(z, x@bestMatch[[i]]@aligns[[j]]@TF@name)})))
					}
					
					if (isInTFNames)
					{
						TFname <- TRUE
#check TF position
						if (j<=as.integer(f[[filterNumber]]@top[[h]]) )
						{
							topX <- TRUE
#check evalue
							if (  x@bestMatch[[i]]@aligns[[j]]@evalue<=as.numeric(f[[filterNumber]]@evalueMax[[h]]) )
							{
								eval <- TRUE
							}
#check motif length
							if( dim(x@bestMatch[[i]]@aligns[[j]]@TF@pwm)[2]>=as.integer(f[[filterNumber]]@lengthMax[[h]]))
							{
								len <- FALSE
							}
							if(TFname && eval && topX && len && name)
							{validMotifs.tmp=c(validMotifs.tmp, j)}
						}
					}
				}
				if (TFname && eval && topX && len && name)
				{
					select[[h]] <- c(select[[h]], i)
					if (is.null(f[[filterNumber]]@valid[[h]]))
					{
						validMotifs <- c(validMotifs, min(validMotifs.tmp))
					}
					else
					{
						validMotifs <- c(validMotifs, f[[filterNumber]]@valid[[h]][countvalidMotifs])
						countvalidMotifs=countvalidMotifs+1
					}
				}
			}
		}
		
		if (length(f[[filterNumber]]@name) > 1)
		{
			selectedMotifs <- c(selectedMotifs, unlist(select)[duplicated(unlist(select)) & unlist(select)!=-1])
		} 
		else 
		{
			selectedMotifs <- c(selectedMotifs, select[[1]][-1])
		}
	}
#create a new motiv object
	for (l in unique(selectedMotifs))
	{
		selectedMatch <- c(selectedMatch, new("matches", name=x@bestMatch[[l]]@name, aligns=x@bestMatch[[l]]@aligns, similarity=x@bestMatch[[l]]@similarity, valid=validMotifs[countSelectedMotifs])	)  	
		countSelectedMotifs=countSelectedMotifs+1
	}
	
	if (!is.null(selectedMatch))
	{
		res <- new("motiv", input=x@input[unique(selectedMotifs)], bestMatch=selectedMatch, argv=x@argv)
		if (verbose)
		{
			message("motivFilter selected ", length(selectedMatch), " motifs.")
		}
	}
	else
	{
		stop ("\t No match.")
	}
	return(res)
}


#####COMBINE#####

combineMotiv <- function (x, filt, name, verbose, exact)
{
	length.motif=length(x@bestMatch)	
	originalNames <- names(x)
	motifsToCombine <- NULL
	
	if (verbose)
	{
		message("motiv object contains ", length.motif, " motifs.")
	}
	
	for (filterNumber in 1:length(filt))
	{
		combinedMotifs <- NULL
		selectedMotifs <- NULL
		n=0
		filter <- filt[[filterNumber]]@filters
		length.filters=length(filter)
		
		for (k in 1:length.filters)
		{
			select <- list()
			for (h in 1:length(filter[[k]]@name))
			{
				select[[h]]=-1
				for (i in 1:length.motif)
				{
					mname <- FALSE
					TFname <- FALSE
					eval <- FALSE
					topX <- FALSE
					len <- TRUE
#check motif name
					if (exact && filter[[k]]@name[[h]]!="")
					{
						isInNames <- any(unlist(lapply(filter[[k]]@name[[h]], function(z){z==x@bestMatch[[i]]@name})))
					} 
					else 
					{
						isInNames <- any(unlist(lapply(filter[[k]]@name[[h]], function(z){grepl(z, x@bestMatch[[i]]@name)})))
					}
					if (isInNames)
					{
						mname <- TRUE
					}
					
					for (j in 1:length(x@bestMatch[[i]]@aligns))
					{
#check TF names
						if (exact && filter[[k]]@tfname[[h]]!="")
						{
							isInTFNames <- any(unlist(lapply(filter[[k]]@tfname[[h]], function(z){z==x@bestMatch[[i]]@aligns[[j]]@TF@name})))
						} 
						else 
						{
							isInTFNames <- any(unlist(lapply(filter[[k]]@tfname[[h]], function(z){grepl(z, x@bestMatch[[i]]@aligns[[j]]@TF@name)})))
						}
						
						if (isInTFNames)
						{
							TFname <- TRUE	
#check TF position
							if (j<=as.integer(filter[[k]]@top[[h]]) )
							{
								topX <- TRUE
#check evalue
								if (  x@bestMatch[[i]]@aligns[[j]]@evalue<=as.numeric(filter[[k]]@evalueMax[[h]]) )
								{
									eval <- TRUE
								}
#check motif length
								if( dim(x@bestMatch[[i]]@aligns[[j]]@TF@pwm)[2]>=as.integer(filter[[k]]@lengthMax[[h]]))
								{
									len <- FALSE
								}
							}
						}
					}
					if (TFname && eval && topX && len && mname)
					{
						select[[h]] <- c(select[[h]], i)
					}
				}
			}
			if (length(filter[[k]]@name) > 1)
			{
				selectedMotifs <- c(selectedMotifs, unlist(select)[duplicated(unlist(select)) & unlist(select)!=-1])
			} 
			else
			{
				selectedMotifs <- c(selectedMotifs, select[[1]][-1])
			}
		}
		if (length(selectedMotifs) >0)
		{
			motifsToCombine <- c(motifsToCombine, selectedMotifs)
			for (l in selectedMotifs)
			{
				n=n+1
				if (!is.null(name))
				{
					combinedMotifs <- c(combinedMotifs, originalNames[l])
					x@bestMatch[[l]]@similarity <- name[filterNumber]
				}
				else
				{
					combinedMotifs <- c(combinedMotifs, originalNames[l])
					x@bestMatch[[l]]@similarity <- paste("s", filterNumber, sep="")
				}
			}
			if(verbose)
			{message(n," motifs combined : ", paste(combinedMotifs, collapse=" "))}
			
		} 
		else
		{
			if(verbose)
			{warning(paste("No motid combined for filter", k,".", sep=""))}
		}
	}
	if(length(motifsToCombine[duplicated(motifsToCombine)==TRUE])>=1)
	{
		warning("Some motifs have been combined many times.")
	}
	return(x)
}

######SPLIT#####

splitMotiv <- function(x, f, drop, exact, verbose)
{
	split <- list()
	motivDrop <- x
	bestMatch <- NULL
	
	if (verbose)
	{message("motiv object contains ", length(x@bestMatch)," motifs.")}
	
	if(length(f)==1)
	{f=as.list(f)}
	
	for (i in 1:length(f))
	{
#get motiv
		split[[i]] <- filter(x, f[[i]], verbose=FALSE, exact=exact)
		names(split)[i] <- paste("filter", i, sep="")			
		motivDrop@input <- motivDrop@input[!names(motivDrop@input)%in%names(split[[i]]@input)]	
		if (verbose)
		{
			message(length(split[[i]]@bestMatch), paste( " motifs selected for 'filter", i, "'", sep=""))
		}
	}
	
	if (length(motivDrop@input)!=0)
	{
#create new motiv objects
		for (l in 1:length(x@bestMatch))
		{
			if (x@bestMatch[[l]]@name %in% unique(names(motivDrop@input)))
			{
				bestMatch <- c(bestMatch, new("matches", name=x@bestMatch[[l]]@name, aligns=x@bestMatch[[l]]@aligns, similarity=x@bestMatch[[l]]@similarity, valid=x@bestMatch[[l]]@valid)	) 
			}
		}
		motivDrop@bestMatch <- bestMatch
		if(!drop)
		{
			split[[i+1]] <- motivDrop
			names(split)[i+1]<-"remaining"
		}
	}
	
	if (verbose )
	{
		if (drop)
		{
			message(length(motivDrop@input), " motifs droped.")
		}
		else
		{
			message( length(motivDrop@input), " remaining motifs selected.")
		}
	}
	return (split)
}
