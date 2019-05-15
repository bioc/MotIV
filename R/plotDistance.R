
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
			if (length(pos1)>=1 & length(pos2)>=1)
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

plotDistanceHistogram<- function ( xlim, pos, distance, column, row, names,  motifName1, motifName2, strands, col, border, nclass, cex)
{
	pushViewport(plotViewport(c(7/length(pos), 0.5,7/length(pos),0.5), name="plot"))
	distance.xlim <- (distance[distance>xlim[1] & distance<xlim[2]])/(2*xlim[2]) +0.5	
	distance.hist=hist(distance.xlim, breaks=seq(0,1,length=nclass+1), plot=F)$count
	do.call("panel.rect", c(list(xright=seq(0,1-1/nclass,length=nclass), ybottom=rep(0,nclass), xleft=seq(1/nclass,1,length=nclass),ytop=distance.hist/(max(distance.hist)+max(distance.hist)/5), col=col[1], border=border)))	
	popViewport() #end plotViewport "plot"
	grid.text(x=unit(0.5,"npc"), y=unit(0.92,"npc"), paste("d(", names[1],"-", names[2],")"), gp=gpar(cex=cex))
}

##############

plotDistanceDensity <- function ( xlim,  pos, distance,  column, row,  names, strands, col, lwd, lty, bw, cex)
{
	pushViewport(plotViewport(c(7/length(pos), 0.5,7/length(pos),0.5), name="plot"))
	distance.xlim <- distance[distance>xlim[1] & distance<xlim[2]]
	density <- do.call ("density", c(list(x=distance.xlim, bw)))
	density.x <- (density$x[density$x>xlim[1] & density$x<xlim[2]] - xlim[1]) /(2*xlim[2]) 
	density.y <- density$y[density$x>xlim[1] & density$x<xlim[2]]
	do.call("panel.lines", c( list(x=density.x, y=density.y/max(density.y), col=col[1], lwd=lwd, lty=lty))) 
	popViewport() #end plotViewport "plot"
	grid.text(x=unit(0.5,"npc"), y=unit(0.92,"npc"), paste("d(", names[1],"-", names[2],")"), gp=gpar(cex=cex))
}
#########

	plotDistanceVenn <- function (cooccurences, motifName1, motifName2, nSequences, cex, vcol)
{
	sequences.n=cooccurences[2, 1]
	sequences1.n=cooccurences[2, 2]
	sequences2.n=cooccurences[1, 1]
	sequences1.alone.n=sequences1.n-sequences.n
	sequences2.alone.n=sequences2.n-sequences.n
	theta <- 2 * pi * (1:360)/360
	grid.lines(unit(0.35+0.3*cos(theta), "npc"), unit(0.5+0.3 * sin(theta), "npc"), gp=gpar(col=vcol[1], lwd=2))
	grid.lines(unit(0.65+0.3*cos(theta), "npc"), unit(0.5+0.3 * sin(theta), "npc"), gp=gpar(col=vcol[2], lwd=2))
	do.call("grid.text", list(x=unit(0.5,"npc"), y=unit(0.9,"npc"), label="# of sequences containing", just="center", gp=gpar(col="black", cex=5*cex/7)))
	do.call("grid.text", list(x=unit(0.1,"npc"), y=unit(0.1,"npc"), label=motifName1, just="left", gp=gpar(col=vcol[1], cex=cex)))
	do.call("grid.text", list(x=unit(0.9,"npc"), y=unit(0.1,"npc"), label=motifName2, just="right", gp=gpar(col=vcol[2], cex=cex)))
	do.call("grid.text", list(x=unit(c(0.5, 0.5),"npc"), y=unit(c(0.52, 0.37),"npc"), label=c(sequences.n , paste(round(100*sequences.n/nSequences,1),"%", sep="")), just="center", gp=gpar(col="black", cex=c(1,0.7)*cex)) )
	do.call("grid.text", list(x=unit(c(0.22,0.24),"npc"), y=unit(c(0.52,0.37),"npc"), label=c(sequences1.alone.n, paste(round(100*sequences1.alone.n/nSequences,1),"%", sep="")), just="center", gp=gpar(col=vcol[1], cex=c(1,0.7)*cex)))
	do.call("grid.text", list(x=unit(c(0.78,0.77),"npc"), y=unit(c(0.52,0.37),"npc"), c(sequences2.alone.n, paste(round(100*sequences2.alone.n/nSequences,1),"%", sep="")), just="center", gp=gpar(col=vcol[2], cex=c(1,0.7)*cex))) 
}

#########

plotDistance <- function (pos, table, strand, main, group, bysim, xlim, nSequences, col, border, lwd, lty, nclass, bw, cex, vcol)
{
	nmotifs=length(pos)
	grid.newpage()
	grid.text(main, gp=gpar(col="black", font=2, cex=1.5), y=unit(0.98,"npc"))
	grid.text("distance", gp=gpar(col="black", font=3, cex=1), y=unit(0.015,"npc"))
	
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
				maxL <- max(mcols(pos[[row]]@positionVector)$lengthPeak, mcols(pos[[distance.col ]]@positionVector)$lengthPeak)

				 if (is.null(xlim))
				 {
					xl <- c(-maxL, maxL)
				 } else {
					xl <- xlim
				 }
		
				distance <- motivDistance(pos[[distance.col]]@positionVector, pos[[row]]@positionVector, dist.commonSeq, motifNames[distance.col], motifNames[row],  c("+","-"), maxL)
				if (!is.null(distance[[1]]))
				{
					plotDistanceDensity(xl, pos, distance[[1]], distance.col, row, distance[[2]], c("+","-"), col, lwd, lty, bw, cex)
					plotDistanceHistogram(xl, pos, distance[[1]], distance.col, row,  distance[[2]], motifNames[distance.col], motifNames[row], c("+","-"), col, border, nclass, cex)
				}
							
				pushViewport(plotViewport(c(0, 0.5,0,0.5), name="axis"))
				panel.segments(0,0.13,1,0.13, col="black", lwd=1)
				panel.segments(x0=c(0.5, 0.25, 0.75), y0=c(0.1, 0.1, 0.1), x1=c(0.5, 0.25, 0.75), y1=c(0.155, 0.155, 0.155), col="black", lwd=1)
				panel.text(x=c(0.25,0.5,0.75), y=0.045, labels=round(c(-xl[2]/2,0, xl[2]/2),0), cex=max(2.5/nmotifs,0.5)*cex)				
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
			grid.text(x=unit(0.5,"npc"), y=unit(0.5,"npc"), pos[[row]]@similarity, gp=gpar(cex=6/nmotifs*cex, srt=45, col="red", font=2))	
		} 
		else 
		{
			grid.text(x=unit(0.5,"npc"), y=unit(0.5,"npc"), pos[[row]]@motifName, gp=gpar(cex=6/nmotifs*cex, srt=45, col="red", font=2))	
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
				plotDistanceVenn (cooccurences, motifNames[row], motifNames[venn.col], nSequences, cex, vcol)
				popViewport() #end vpseqcom	
			}
		}
	}
	popViewport() #end vp
	popViewport() #end plotViewport "grid"
}
