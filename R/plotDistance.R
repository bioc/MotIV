
#####DISTANCE#####
motivDistance <- function (positionVector1, positionVector2, commonSequences, name1, name2, strands, sequencesLength)
{
positions1=as.data.frame(positionVector1)
positions2=as.data.frame(positionVector2)
distance <- NULL
#Apply the method depending of the lowest standard deviation
 if(var(positions1$peakPos) <= var(positions2$peakPos) )
    {
		for (seqID in commonSequences)
		{
			pos1=unique(positions1$peakPos[  positions1$seqID==seqID & positions1$strand%in%strands ])
			pos2=unique(positions2$peakPos[ positions2$seqID==seqID & positions2$strand%in%strands ])
			if (length(pos1)>=1& length(pos2)>=1)
			{
				pos1.center=pos1[which(abs(pos1-sequencesLength)==min(abs(pos1-sequencesLength)))] #select center motif
				distance <- c(distance, pos1.center-pos2) #distance
			}
		}
		names <- c(name1, name2)
	}	 else
    {
		for (seqID in commonSequences)
		{
			pos1=unique(positions2$peakPos[ positions2$seqID==seqID & positions2$strand%in%strands ])
			pos2=unique(positions1$peakPos[  positions1$seqID==seqID & positions1$strand%in%strands ])
			if (length(pos1)>=1& length(pos2)>=1)
			{
				pos1.center=pos1[which(abs(pos1-sequencesLength)==min(abs(pos1-sequencesLength)))] #select center motif
				distance <- c(distance, pos1.center-pos2) #distance
			}
		}
		names <- c(name2, name1)
    }
	return (list(distance, names))
}


#########

plotDistanceHistogram<- function (harg, xlim, strand, pos, distance, column, row, names,  motifName1, motifName2, strands, sequencesLength)
{
	xlength=xlim[2]
	nint = harg[["nint"]]
	pushViewport(plotViewport(c(7/length(pos), 0.5,7/length(pos),0.5), name="plot"))
		
	# if (!strand)
	# {
	#	distance <- motivDistance(pos[[column]]@positionVector, pos[[row]]@positionVector, commonSequences,  motifName1, motifName2, strands, sequencesLength)
		distance.xlim <- (distance[distance>xlim[1] & distance<xlim[2]])/(2*xlength) +0.5	
		distance.hist=hist(distance.xlim, breaks=seq(0,1,length=nint+1), plot=F)$count
		do.call("panel.rect", c(list(xright=seq(0,1-1/nint,length=nint), ybottom=rep(0,nint), xleft=seq(1/nint,1,length=nint),ytop=distance.hist/(max(distance.hist)+max(distance.hist)/5)), harg))
	# } else {
		# distance<- list()
		# distance.hist <- list()
		# for (i in 1:2)	
		# {
			# distance <- motivDistance(pos[[column]]@positionVector, pos[[row]]@positionVector, commonSequences,  motifName1, motifName2, strands[i], sequencesLength)
			# distance.xlim <- (distance[[1]][distance[[1]]>xlim[1] & distance[[1]]<xlim[2]])/xlength
			# distance.hist[[i]]=hist(distance.xlim+0.5, breaks=seq(0,1,length=nint+1), plot=F)$count
		# }		
		# arg <- harg[names(harg)!="border"] 
		# hist.max <- max(distance.hist[[1]], distance.hist[[2]])
		# do.call("panel.rect", c(list(xright=seq(0,1-1/nint,length=nint), ybottom=rep(0,nint), xleft=seq(1/nint,1,length=nint),ytop=distance.hist[[1]]/(hist.max+hist.max/5), border=2), arg))
		# do.call("panel.rect", c(list(xright=seq(0,1-1/nint,length=nint), ybottom=rep(0,nint), xleft=seq(1/nint,1,length=nint),ytop=distance.hist[[2]]/(hist.max+hist.max/5), border=4), arg))
	# }
	popViewport() #end plotViewport "plot"
	grid.text(x=unit(0.5,"npc"), y=unit(0.92,"npc"), paste("d(", names[1],"-", names[2],")"))
}

##############

plotDistanceDensity <- function (darg, carg, xlim, strand, pos, distance,  column, row,  names, motifName1, motifName2, strands, sequencesLength)
{
	xlength=xlim[2]
	pushViewport(plotViewport(c(7/length(pos), 0.5,7/length(pos),0.5), name="plot"))
	# if(!strand)
	# {		
		#distance <- motivDistance(pos[[column]]@positionVector, pos[[row]]@positionVector, commonSequences, motifName1, motifName2, c("+","-"), sequencesLength)
		distance.xlim <- distance[distance>xlim[1] & distance<xlim[2]]
		density <- do.call ("density", c(list(x=distance.xlim),darg))
		density.x <- (density$x[density$x>xlim[1] & density$x<xlim[2]] - xlim[1])/(2*xlength)
		density.y <- density$y[density$x>xlim[1] & density$x<xlim[2]]
		do.call("panel.lines", c( list(x=density.x, y=density.y/max(density.y)), carg)) 
	# } else {
		# density.x <- list(NULL, NULL)
		# density.y <- list(NULL, NULL)
		# for (i in 1:2)
		# {
			# distance <- motivDistance(pos[[column]]@positionVector, pos[[row]]@positionVector, commonSequences, motifName1, motifName2, strands[i], sequencesLength)
			# distance.xlim <- distance[[1]][distance[[1]]>xlim[1] & distance[[1]]<xlim[2]]
			# if (length(distance)>1)
			# {
				# density <- do.call ("density", c(list(x=distance.xlim),darg))
				# density.x[[i]] <- (density$x[density$x>xlim[1] & density$x<xlim[2]] -xlim[1])/xlength
				# density.y[[i]] <- density$y[density$x>xlim[1] & density$x<xlim[2]]
			# }
		# }
		# if (length(density.y[[1]])>1 )
		# {
			# carg[["col"]] <- "red"
			# do.call("panel.lines", c( list(x=density.x[[1]], y=density.y[[1]]/max(density.y[[1]], density.y[[2]])), carg)) 
		# }
		# if (length(density.y[[2]])>1 )
		# {
			# carg[["col"]] <- "blue"
			# do.call("panel.lines", c( list(x=density.x[[2]], y=density.y[[2]]/max(density.y[[1]], density.y[[2]])), carg)) 
		# }		
	# }
	popViewport() #end plotViewport "plot"
	grid.text(x=unit(0.5,"npc"), y=unit(0.92,"npc"), paste("d(", names[1],"-", names[2],")"))
}
#########

plotDistanceVenn <- function (cooccurences,  varg, motifName1, motifName2, nSequences)
{
	sequences.n=cooccurences[2, 1]
	sequences1.n=cooccurences[2, 2]
	sequences2.n=cooccurences[1, 1]
	sequences1.alone.n=sequences1.n-sequences.n
	sequences2.alone.n=sequences2.n-sequences.n
	do.call("grid.circle", list(x=unit(0.4,"npc"), y=unit(0.5,"npc"), r=unit(0.3,"npc"), gp=gpar(col=varg[["col"]][1], lwd=varg[["lwd"]])))
	do.call("grid.circle", list(x=unit(0.6,"npc"), y=unit(0.5,"npc"), r=unit(0.3,"npc"), gp=gpar(col=varg[["col"]][2], lwd=varg[["lwd"]])))
	do.call("grid.text", list(x=unit(0.1,"npc"), y=unit(0.9,"npc"), label=motifName1, just="left", gp=gpar(col=varg[["col"]][1], cex=varg[["cex"]])))
	do.call("grid.text", list(x=unit(0.9,"npc"), y=unit(0.9,"npc"), label=motifName2, just="right", gp=gpar(col=varg[["col"]][2], cex=varg[["cex"]])))
	do.call("grid.text", list(x=unit(c(0.5, 0.5),"npc"), y=unit(c(0.5, 0.60),"npc"), label=c(sequences.n , paste(round(100*sequences.n/nSequences,1),"%", sep="")), just="center", gp=gpar(col="black", cex=c(1,0.7)*varg[["cex"]])) )
	do.call("grid.text", list(x=unit(c(0.35,0.20),"npc"), y=unit(c(0.5,0.75),"npc"), label=c(sequences1.alone.n, paste(round(100*sequences1.alone.n/nSequences,1),"%", sep="")), just="center", gp=gpar(col=varg[["col"]][1], cex=c(1,0.7)*varg[["cex"]])))
	do.call("grid.text", list(x=unit(c(0.65,0.80),"npc"), y=unit(c(0.5,0.75),"npc"), c(sequences2.alone.n, paste(round(100*sequences2.alone.n/nSequences,1),"%", sep="")), just="center", gp=gpar(col=varg[["col"]][2], cex=c(1,0.7)*varg[["cex"]]))) 
	
	#do.call("grid.text", list(x=unit(c(0.35,0.35),"npc"), y=unit(c(0.5,0.35),"npc"), label=c(sequences1.alone.n, paste(round(100*sequences1.alone.n/sequences1.n,1),"%", sep="")), just="center", gp=gpar(col=varg[["col"]][1], cex=c(1,0.7)*varg[["cex"]])))
	#do.call("grid.text", list(x=unit(c(0.65,0.65),"npc"), y=unit(c(0.5,0.35),"npc"), c(sequences2.alone.n, paste(round(100*sequences2.alone.n/sequences2.n,1),"%", sep="")), just="center", gp=gpar(col=varg[["col"]][2], cex=c(1,0.7)*varg[["cex"]])))   
}

#########

plotDistance <- function (pos, table, strand, main, group, bysim, xlim, sequencesLength, nSequences, harg, darg, carg, varg,...)
{
	nmotifs=length(pos)
	grid.newpage()
	grid.text(main, gp=gpar(col="black", font=2, cex=1.5), y=unit(0.98,"npc"))
	grid.text("distance", gp=gpar(col="black", font=3, cex=1), y=unit(0.015,"npc"))
	if (strand)
	{
		grid.text("Forward", gp=gpar( font=3, cex=0.7,  col="blue"), y=unit(0.02,"npc"), x=unit(0.05,"npc"), just="top")
		grid.text("Reverse", gp=gpar( font=3, cex=0.7,  col="red"), y=unit(0.02,"npc"), x=unit(0.15,"npc"), just="top")
	}
	if(bysim)
	{
		motifNames <- similarity(pos)
	} else {
		motifNames <- sapply(pos, function(x)x@motifName)
	}
	
	pushViewport(plotViewport(c(7/nmotifs, 0,7/nmotifs,0), name="grid"))
	grid.segments(c(1:(nmotifs-1)/nmotifs, rep(0, nmotifs-1)), c(rep(0, nmotifs-1), 1:(nmotifs-1)/nmotifs), c(1:(nmotifs-1)/nmotifs, rep(1, nmotifs-1)), c(rep(1, nmotifs-1), 1:(nmotifs-1)/nmotifs), gp = gpar(col = "grey90", lwd=1, lty=2))
	grid.segments(c(1:(nmotifs-1)/nmotifs, rep(0, nmotifs-1)), c((nmotifs-1):1/nmotifs, 1:(nmotifs-1)/nmotifs), c(1:(nmotifs-1)/nmotifs, (nmotifs-1):1/nmotifs), c(rep(1, nmotifs-1), 1:(nmotifs-1)/nmotifs), gp = gpar(col = "grey90", lwd=2))
	
	vp <- viewport(layout=grid.layout(nmotifs, nmotifs)) #parse case
	pushViewport(vp)
	par(mar=c(2,1,1,1))
	if (is.null(xlim))
	{xlim <- c(-sequencesLength, sequencesLength)}
	xlength=xlim[2]
	for (row in 1:nmotifs)
	{
#plot distance
		if (row>1)
		{
			for (distance.col in 1:(row-1))
			{	
				vpdistance <- viewport(layout.pos.col=distance.col, layout.pos.row=row) #seqLogo
				pushViewport(vpdistance)	
				dist.commonSeq <- which(table[,distance.col]!=0 & table[,row]!=0)

				
				distance <- motivDistance(pos[[distance.col]]@positionVector, pos[[row]]@positionVector, dist.commonSeq, motifNames[distance.col], motifNames[row],  c("+","-"), sequencesLength)

				if(is.null(carg[["plot"]]) || carg[["plot"]])
				{ plotDistanceDensity( darg, carg, xlim, strand, pos, distance[[1]], distance.col, row, distance[[2]], motifNames[distance.col], motifNames[row], c("+","-"), sequencesLength)}
				if (is.null(harg[["plot"]]) || harg[["plot"]])
				{plotDistanceHistogram(harg, xlim, strand, pos, distance[[1]], distance.col, row,  distance[[2]], motifNames[distance.col], motifNames[row], c("+","-"), sequencesLength)}
				
				pushViewport(plotViewport(c(0, 0.5,0,0.5), name="axis"))
				panel.segments(0,0.13,1,0.13, col="black", lwd=1)
				panel.segments(x0=c(0.5, 0.25, 0.75), y0=c(0.1, 0.1, 0.1), x1=c(0.5, 0.25, 0.75), y1=c(0.155, 0.155, 0.155), col="black", lwd=1)
				panel.text(x=c(0.25,0.5,0.75), y=0.045, labels=round(c(-xlength/2,0, xlength/2),0), cex=2.5/nmotifs)				
				# panel.text(x=0.5, y=0.045, labels="0", cex=2.5/nmotifs)
				# panel.text(x=0.25, y=0.045, labels=round(*xlength/2), cex=2.5/nmotifs)
				# panel.text(x=0.75, y=0.045, labels=round(xlength/2), cex=2.5/nmotifs)
				popViewport() #end plotViewport "axis"
				popViewport() #end vpdistance
			}
		}
		
#print motif name
		vpname <- viewport(layout.pos.col=row, layout.pos.row=row) #seqLogo
		pushViewport(vpname)	
		pushViewport(	viewport(x=0.5, y=0.5, width=0.6, height=0.25, angle=45, name="rect"))
		grid.rect(x=unit(0.5,"npc"), y=unit(0.5,"npc"), height=unit(1,"npc"), width=unit(1,"npc"), gp=gpar(col="grey95", fill="grey95"))				
		
		if (bysim)
		{
			grid.text(x=unit(0.5,"npc"), y=unit(0.5,"npc"), pos[[row]]@similarity, gp=gpar(cex=6/nmotifs, srt=45, col="red", font=2))	
		} 
		else 
		{
			grid.text(x=unit(0.5,"npc"), y=unit(0.5,"npc"), pos[[row]]@motifName, gp=gpar(cex=6/nmotifs, srt=45, col="red", font=2))	
		}
		popViewport()	 #end "rect"
		popViewport()	#end vpname
		
#plot common sequences
		if (row<nmotifs)
		{
			for(venn.col in (row+1):nmotifs)
			{
				vpseqcom <- viewport(layout.pos.col=venn.col, layout.pos.row=row) #seqLogo
				pushViewport(vpseqcom)	
				venn.commonSeq <- which(table[,venn.col]!=0 & table[,row]!=0)
				cooccurences = cooccurences(table[,c(venn.col, row)])
				plotDistanceVenn (cooccurences, varg, motifNames[row], motifNames[venn.col], nSequences)
				popViewport() #end vpseqcom	
			}
		}
	}
	popViewport() #end vp
	popViewport() #end plotViewport "grid"
}
