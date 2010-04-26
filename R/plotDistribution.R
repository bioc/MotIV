

plotDistributionHistogram <- function (pos, curr, size,  strand, harg, sequencesLength)
{
	nint = harg[["nint"]]
	if (!strand)
	{
	#	pos.mean <- (end(pos[[curr]]@positionVector) + start(pos[[curr]]@positionVector))/(2*sequencesLength)
	pos.mean <- pos[[curr]]@positionVector$peakPos/sequencesLength
		distance.hist=hist(pos.mean[pos.mean>0 & pos.mean<1], breaks=seq(0,1,length=nint+1), plot=F)$count
		do.call("panel.rect", c(list(xright=seq(0,1-1/nint,length=nint), ybottom=rep(0,nint), xleft=seq(1/nint,1,length=nint),ytop=distance.hist/(max(distance.hist)+max(distance.hist)/5)), harg))
	} else {
		strands = c("+", "-")
		dist<- list()
		distance.hist <- list()
		for (i in 1:2)
		{
		#	pos.mean <- (end(pos[[curr]]@positionVector)[pos[[curr]]@positionVector$gadem.strand==strands[i]] + start(pos[[curr]]@positionVector)[pos[[curr]]@positionVector$gadem.strand==strands[i]])/(2*sequencesLength)
		pos.mean <- (pos[[curr]]@positionVector$peakPos[pos[[curr]]@positionVector$strand==strands[i]] + pos[[curr]]@positionVector$peakPos[pos[[curr]]@positionVector$strand==strands[i]])/(2*sequencesLength)
			distance.hist[[i]]=hist(pos.mean , breaks=seq(0,1,length=nint+1), plot=F)$count
		}		
		arg <- harg[names(harg)!="border"]
		hist.max <- max(distance.hist[[1]], distance.hist[[2]])
		do.call("panel.rect", c(list(xright=seq(0,1-1/nint,length=nint), ybottom=rep(0,nint), xleft=seq(1/nint,1,length=nint),ytop=distance.hist[[1]]/(hist.max+hist.max/5), border=2), arg))
		do.call("panel.rect", c(list(xright=seq(0,1-1/nint,length=nint), ybottom=rep(0,nint), xleft=seq(1/nint,1,length=nint),ytop=distance.hist[[2]]/(hist.max+hist.max/5), border=4), arg))
	}
}

#######

plotDistributionDensity <- function(pos, curr, size,  strand, darg, carg, sequencesLength)
{
	xlim=c(0,sequencesLength)
	if(!strand)
	{
		#pos.mean <- (end(pos[[curr]]@positionVector) + start(pos[[curr]]@positionVector))/2
		pos.mean  <- pos[[curr]]@positionVector$peakPos
		density <- do.call ("density", c(list(x=pos[[curr]]@positionVector$peakPos),darg))
		density.x <- (density$x[density$x>xlim[1] & density$x<xlim[2]])/sequencesLength
		density.y <- density$y[density$x>xlim[1] & density$x<xlim[2]]
		do.call("panel.lines", c(list(x=density.x, y=density.y/max(density.y)), carg))
#pos.mean.center <- pos.mean-mean(pos.mean)
# kurtosis=mean(pos.mean.center ^4)/mean(pos.mean.center ^2)^2-3
# grid.text(paste("Kurtosis\n", round(kurtosis,3)), x=unit(-0.2,"npc"), y=unit(0.6,"npc"), gp=gpar(cex=max(0.6,2/size), font=2))
		grid.text(paste("sd\n", round(sd(pos.mean),3)), x=unit(1.2,"npc"), y=unit(0.6,"npc"), gp=gpar(cex=max(0.6,2/size), font=2))
	} else {
		strands = c("+", "-")
		distance.x <- list(NULL, NULL)
		distance.y <- list(NULL, NULL)
		for (i in 1:2)
		{
			pos.mean <- (end(pos[[curr]]@positionVector)[pos[[curr]]@positionVector$strand==strands[i]] + start(pos[[curr]]@positionVector)[pos[[curr]]@positionVector$strand==strands[i]])/2
			if (length(pos.mean)>1)
			{
				density <- do.call ("density", c(list(x=pos.mean),darg))
				distance.x[[i]] <- (density$x[density$x>xlim[1] & density$x<xlim[2]] -xlim[1])/sequencesLength
				distance.y[[i]] <- density$y[density$x>xlim[1] & density$x<xlim[2]]
			}
		}
		if (length(distance.y[[1]])>1 )
		{
			carg[["col"]] <- "red"
			do.call("panel.lines", c( list(x=distance.x[[1]], y=distance.y[[1]]/max(distance.y[[1]], distance.y[[2]])), carg)) 
		}
		if (length(distance.y[[2]])>1 )
		{
			carg[["col"]] <- "blue"
			do.call("panel.lines", c( list(x=distance.x[[2]], y=distance.y[[2]]/max(distance.y[[1]], distance.y[[2]])), carg)) 
		}	
	}
}




#####DISTRIBUTION#####
plotDistribution <- function (pos, group, main, sort, ncol, nrow, strand, bysim, sequencesLength, harg, darg, carg)
{
	positionVector <- list()
	for (i in 1:length(pos))
	{
		positionVector[[i]]=(end(pos[[i]]@positionVector)+ start(pos[[i]]@positionVector))/2
	}
	
	var <- sapply(positionVector, function(x){var(as.integer(x))})
	
	if (group) 
	{
		names(var)=unique(similarity(pos))
	} 
	else 
	{
		names(var)=sapply(pos, function(x){x@motifName})
	}
	
	if(sort)
	{
		var.sort=sort(var)
	} 
	else 
	{
		var.sort=var
	}
	if(ncol==0 && nrow==0)
	{
		ncol=ceiling(sqrt(length(pos)))
		nrow=ncol
	} 
	else	if (ncol==0 && nrow!=0)
	{
		ncol=ceiling(length(pos)/nrow)
	} 
	else if ((ncol!=0 && nrow==0))
	{
		nrow=ceiling(length(pos)/ncol)
	}
	
	if (nrow*ncol< length(pos) )
	{
		print(paste("Warning : layout is not sufficient, some results couldn't be plot."))
	}
	layout <- c(nrow, ncol)
	size=max(layout[1], layout[2])	
	
	grid.newpage()
	grid.text(main, gp=gpar(col="black", font=2, cex=1.5), y=unit(0.98,"npc"))
	grid.text("position", gp=gpar(col="black", font=3, cex=1), y=unit(0.01,"npc"))
	if (strand)
	{
		grid.text("Forward", gp=gpar( font=3, cex=0.7,  col="blue"), y=unit(0.02,"npc"), x=unit(0.18,"npc"), just="top")
		grid.text("Reverse", gp=gpar( font=3, cex=0.7,  col="red"), y=unit(0.02,"npc"), x=unit(0.27,"npc"), just="top")
	}
	pushViewport(plotViewport(c(7/layout[1], 0,7/layout[1],0), name="grid"))
	vp  <- viewport(layout=grid.layout(layout[1], layout[2]))	#1rst lvl
	pushViewport(vp)
	grid.segments(c(1:(layout[2]-1)/layout[2], rep(0, layout[1]-1)), c(rep(0, layout[2]-1), 1:(layout[1]-1)/layout[1]), c(1:(layout[2]-1)/layout[2], rep(1, layout[1]-1)), c(rep(1, layout[2]-1), 1:(layout[1]-1)/layout[1]), gp = gpar(col = "grey90", lwd=2))
	
#plot Logo
	countMotifs=1
	for (row in 1:layout[1])
	{
		for (column in 1:layout[2])
		{
			if (countMotifs<=length(pos))
			{
				curr=which(names(var.sort)[countMotifs]==names(var))
				
				vpcase <- viewport(layout.pos.col=column, layout.pos.row=row)
				pushViewport(vpcase)
				
				plotMotif ( pos, curr, column,  3, layout, bysim, group)
				
				grid.text("Distribution", x=unit(0,"npc"), y=unit(-0.15,"npc"), gp=gpar(cex=max(0.6,2.5/size), font=2))
				popViewport() #end vpmotifrev

#plot distribution
				vpdistribution <- viewport(layout.pos.col=2:3, layout.pos.row=2:4) #seqLogo
				pushViewport(vpdistribution)		
				pushViewport(plotViewport(c(max(5/size,0.9), 0,3.5/size,0), name="plot"))
				
				if(is.null(carg[["plot"]]) || carg[["plot"]])
				{plotDistributionDensity(pos, curr, size, strand, darg, carg, sequencesLength)}
				if (is.null(harg[["plot"]]) || harg[["plot"]])
				{plotDistributionHistogram(pos, curr, size,  strand, harg, sequencesLength)}
				
				popViewport() #end vpdistribution
				panel.segments(0,0.13,1,0.13, col="black", lwd=1)
				panel.segments(x0=c(0.25,0.5,0.75), y0=0.1+c(0,0,0), x1=c(0.25,0.5,0.75), y1=0.115+c(0.05,0.05,0.05), col="black", lwd=1)	
				panel.text(x=c(0.25,0.5,0.75), y=0.045, labels=round(c(-sequencesLength/2,0, sequencesLength/2),0), cex=2/size)
				popViewport() #end  plotViewport "plot"
				popViewport() #end vpin
				popViewport() #end vpcase
				popViewport() #end plotViewport "case"
				countMotifs=countMotifs+1
			}
		}
	}
	popViewport() #end vp
	popViewport() #end plotViewport "grid"
}
