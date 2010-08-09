

plotDistributionHistogram <- function (pos, curr, size, strand, maxL, col, border, nclass)
{
	if (!strand)
	{
		pos.mean <- pos[[curr]]@positionVector$peakPos/maxL + 0.5
		distance.hist=hist(pos.mean[pos.mean>0 & pos.mean<1], breaks=seq(0,1,length=nclass+1), plot=F)$count
		do.call("panel.rect", c(list(xright=seq(0,1-1/nclass,length=nclass), ybottom=rep(0,nclass), xleft=seq(1/nclass,1,length=nclass),ytop=distance.hist/(max(distance.hist)+max(distance.hist)/5), col=col[1], border=border)))
	} else {
		strands = c("+", "-")
		dist<- list()
		distance.hist <- list()
		for (i in 1:2)
		{
			pos.mean <- (pos[[curr]]@positionVector$peakPos[pos[[curr]]@positionVector$strand==strands[i]] )/(maxL ) + 0.5
			distance.hist[[i]]=hist(pos.mean , breaks=seq(0,1,length=nclass+1), plot=F)$count
		}		
		hist.max <- max(distance.hist[[1]], distance.hist[[2]])
		do.call("panel.rect", c(list(xright=seq(0,1-1/nclass,length=nclass), ybottom=rep(0,nclass), xleft=seq(1/nclass,1,length=nclass),ytop=distance.hist[[1]]/(hist.max+hist.max/5), border=border[1],col=col[1])))
		do.call("panel.rect", c(list(xright=seq(0,1-1/nclass,length=nclass), ybottom=rep(0,nclass), xleft=seq(1/nclass,1,length=nclass),ytop=distance.hist[[2]]/(hist.max+hist.max/5), border=border[2], col=col[2])))
	}
}

#######

plotDistributionDensity <- function(pos, curr, size,  strand, maxL, col, lwd, lty, bw, cex)
{
	xlim=c(0,maxL)
	if(!strand)
	{
		pos.mean <- pos[[curr]]@positionVector$peakPos + maxL/2 
		density <- do.call ("density", c(list(x=pos.mean),bw=bw))
		density.x <- (density$x[density$x>xlim[1] & density$x<xlim[2]]) /maxL
		density.y <- density$y[density$x>xlim[1] & density$x<xlim[2]]
		do.call("panel.lines", c(list(x=density.x, y=density.y/max(density.y), col=col[1], lw=lwd, lty=lty)))
		grid.text(paste("sd\n", round(sd(pos.mean),3)), x=unit(1.2,"npc"), y=unit(0.6,"npc"), gp=gpar(cex=max(0.6,2/size)*cex, font=2))
	} else {
		strands = c("+", "-")
		distance.x <- list(NULL, NULL)
		distance.y <- list(NULL, NULL)
		for (i in 1:2)
		{
			pos.mean <- (pos[[curr]]@positionVector$peakPos[pos[[curr]]@positionVector$strand==strands[i]] ) +  maxL/2 
			
			if (length(pos.mean)>1)
			{
				density <- do.call ("density", c(list(x=pos.mean, bw=bw)))
				distance.x[[i]] <- (density$x[density$x>xlim[1] & density$x<xlim[2]] -xlim[1])/maxL
				distance.y[[i]] <- density$y[density$x>xlim[1] & density$x<xlim[2]]
			}
		}
		if (length(distance.y[[1]])>1 )
		{
			do.call("panel.lines", c( list(x=distance.x[[1]], y=distance.y[[1]]/max(distance.y[[1]], distance.y[[2]]), col=col[1], lwd=lwd, lty=lty))) 
		}
		if (length(distance.y[[2]])>1 )
		{
			do.call("panel.lines", c( list(x=distance.x[[2]], y=distance.y[[2]]/max(distance.y[[1]], distance.y[[2]]), col=col[2], lw=lwd, lty=lty))) 
		}	
	}
}

#####DISTRIBUTION#####
plotDistribution <- function (pos, group, main, sort, ncol, nrow, strand, bysim, trim, col, border, lwd, lty, nclass, bw, cex)
{
	positionVector <- list()
	for (i in 1:length(pos))
	{
		positionVector[[i]]=pos[[i]]@positionVector[["peakPos"]] 
	}
	
	var <- sapply(positionVector, function(x){var(as.integer(x))})
	
	if (group) 
	{
		names(var)=unique(similarity(pos))
	} 	else 	{
		names(var)=sapply(pos, function(x){x@motifName})
	}
	
	if(sort)
	{
		var.sort=sort(var)
	} 	else 	{
		var.sort=var
	}
	if(ncol==0 && nrow==0)
	{
		ncol=ceiling(sqrt(length(pos)))
		nrow=ncol
	} 	else	if (ncol==0 && nrow!=0)
	{
		ncol=ceiling(length(pos)/nrow)
	} 	else if ((ncol!=0 && nrow==0))
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
	#grid.text(main, gp=gpar(col="black", font=2, cex=1.5), y=unit(0.98,"npc"))
	grid.text("position", gp=gpar(col="black", font=3, cex=cex), y=unit(0.01,"npc"))
	if (strand)
	{
		grid.text("Forward", gp=gpar( font=3, cex=0.7*cex,  col=col[1]), y=unit(0.02,"npc"), x=unit(0.18,"npc"), just="top")
		grid.text("Reverse", gp=gpar( font=3, cex=0.7*cex,  col=col[2]), y=unit(0.02,"npc"), x=unit(0.27,"npc"), just="top")
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
				plotMotif ( pos, curr, column,  3, layout, bysim, group, trim, cex=cex)				
				grid.text("Distribution", x=unit(0,"npc"), y=unit(-0.15,"npc"), gp=gpar(cex=max(0.6,2.5/size)*cex, font=2))
				popViewport() #end vpmotifrev
#plot distribution
				vpdistribution <- viewport(layout.pos.col=2:3, layout.pos.row=2:4) #seqLogo
				pushViewport(vpdistribution)		
				pushViewport(plotViewport(c(max(5/size,0.9), 0,3.5/size,0), name="plot"))
				
				maxL <- max(pos[[curr]]@positionVector$lengthPeak)
				plotDistributionDensity(pos, curr, size, strand, maxL, col, lwd, lty, bw, cex)
				plotDistributionHistogram(pos, curr, size,  strand, maxL, col, border, nclass)
												
				popViewport() #end vpdistribution
				panel.segments(0,0.13,1,0.13, col="black", lwd=1)
				panel.segments(x0=c(0.25,0.5,0.75), y0=0.1+c(0,0,0), x1=c(0.25,0.5,0.75), y1=0.115+c(0.05,0.05,0.05), col="black", lwd=1)	
				panel.text(x=c(0.25,0.5,0.75), y=0.045, labels=round(c(-maxL/4,0, maxL/4),0), cex=2*cex/size)
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
	grid.text(main, gp=gpar(col="black", font=2, cex=max(0.8,3/size)*cex), y=unit(0.999,"npc"), just="top")
}
