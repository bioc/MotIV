
#####DISTANCE#####
motivDistance <- function (positionVector1, positionVector2, commonSequences, method, name1, name2, strand, meanLength)
{
  distance <- NULL
  if (method==1) #Every motifs of a sequence are used to compute the distance between motif 1 and 2
  {
    for (seq in commonSequences)
    {
      pos1=unique(positionVector1$start[  positionVector1$seq==seq & positionVector1$gadem.strand%in%strand ])
      pos2=unique(positionVector2$start[ positionVector2$seq==seq & positionVector2$gadem.strand%in%strand ])
      if (length(pos1)>=1& length(pos2)>=1)
      {
        pos1.center=pos1[which(abs(pos1-meanLength)==min(abs(pos1-meanLength)))] #select center motif
        distance <- c(distance, pos1.center-pos2) #distance
      }
    }
    names <- c(name1, name2)
  }
  else if (method==2)
  {		#Same as 1 but it return the distance of the motif 2 versus motif 1
    for (seq in commonSequences)
    {
      pos1=unique(positionVector2$start[ positionVector2$seq==seq & positionVector2$gadem.strand%in%strand ])
      pos2=unique(positionVector1$start[  positionVector1$seq==seq & positionVector1$gadem.strand%in%strand ])
      if (length(pos1)>=1& length(pos2)>=1)
      {
        pos1.center=pos1[which(abs(pos1-meanLength)==min(abs(pos1-meanLength)))] #select center motif
        distance <- c(distance, pos1.center-pos2) #distance
      }
    }
    names <- c(name2, name1)
  } 
  else if (method==3)
  {		#Apply the method 4 or 5 depending of the lowest standard deviation
    if(var(positionVector1$start) <= var(positionVector2$start) )
    {
      for (seq in commonSequences)
      {
        pos1=unique(positionVector1$start[  positionVector1$seq==seq & positionVector1$gadem.strand%in%strand ])
        pos2=unique(positionVector2$start[ positionVector2$seq==seq & positionVector2$gadem.strand%in%strand ])
        if (length(pos1)>=1& length(pos2)>=1)
        {
          pos1.center=pos1[which(abs(pos1-meanLength)==min(abs(pos1-meanLength)))] #select center motif
          distance <- c(distance, pos1.center-pos2) #distance
        }
      }
      names <- c(name1, name2)
    }
    else
    {
      for (seq in commonSequences)
      {
        pos1=unique(positionVector2$start[ positionVector2$seq==seq & positionVector2$gadem.strand%in%strand ])
        pos2=unique(positionVector1$start[  positionVector1$seq==seq & positionVector1$gadem.strand%in%strand ])
        if (length(pos1)>=1& length(pos2)>=1)
        {
          pos1.center=pos1[which(abs(pos1-meanLength)==min(abs(pos1-meanLength)))] #select center motif
          distance <- c(distance, pos1.center-pos2) #distance
        }
      }
      names <- c(name2, name1)
    }

  } 
  else if (method==4) 
  { #get the motif 1 in the middle and computes relative distance
    for (seq in commonSequences)
    {
      pos1=unique(positionVector2$start[ positionVector2$seq==seq & positionVector2$gadem.strand%in%strand ])
      pos2=unique(positionVector1$start[  positionVector1$seq==seq & positionVector1$gadem.strand%in%strand ])
      if (length(pos1)>=1& length(pos2)>=1)
      {distance <- c(distance, rep(pos1 , length(pos2))-rep(pos2, length(pos1)))} #distance
    }	
    names <- c(name2, name1)
  } 
  else if (method==5) 
  { #Same as 4 but with the motif 2 centered.
    for (seq in commonSequences)
    {
      pos1=unique(positionVector1$start[  positionVector1$seq==seq & positionVector1$gadem.strand%in%strand ])
      pos2=unique(positionVector2$start[ positionVector2$seq==seq & positionVector2$gadem.strand%in%strand ])
      if (length(pos1)>=1& length(pos2)>=1)
      {distance <- c(distance, rep(pos1, length(pos2))-rep(pos2, length(pos1)))} #distance
    }	
    names <- c(name1, name2)
  } 
  else if (method==6) 
  { #get the motif 1 in the middle and computes relative distance for every other motifs
    for (seq in commonSequences)
    {
      pos2=unique(positionVector2$start[ positionVector2$seq==seq & positionVector2$gadem.strand%in%strand ])
      pos1=unique(positionVector1$start[  positionVector1$seq==seq & positionVector1$gadem.strand%in%strand ])
      if (length(pos1)>=1& length(pos2)>=1)
      {
        dist= rep(pos1, length(pos2)) - rep (pos2, length(pos1))
        dist.min = which(abs(dist)==min(abs(dist)))
        pos1.center=pos1[which(abs(pos1-meanLength)==min(abs(pos1-meanLength)))] #select center motif
        distance <- c(distance, dist[dist.min]) #distance
      }
    }
    names <- c(name1, name2)
  }
  return (list(distance, names))
}


plotDistance <- function (pos, strand, main, sub, group, bysim, xlim, method, bw, meanLength)
{
  nmotif=length(pos)
  grid.newpage()
  grid.text(main, gp=gpar(col="black", font=2, cex=1.5), y=unit(0.98,"npc"))
  grid.text("distance", gp=gpar(col="black", font=3, cex=1), y=unit(0.015,"npc"))
  if (strand)
  {
    grid.text("Forward", gp=gpar( font=3, cex=0.7,  col="blue"), y=unit(0.02,"npc"), x=unit(0.05,"npc"), just="top")
    grid.text("Backward", gp=gpar( font=3, cex=0.7,  col="red"), y=unit(0.02,"npc"), x=unit(0.15,"npc"), just="top")
  }
  pushViewport(plotViewport(c(7/nmotif, 0,7/nmotif,0), name="grid"))
  grid.segments(c(1:(nmotif-1)/nmotif, rep(0, nmotif-1)), c(rep(0, nmotif-1), 1:(nmotif-1)/nmotif), c(1:(nmotif-1)/nmotif, rep(1, nmotif-1)), c(rep(1, nmotif-1), 1:(nmotif-1)/nmotif), gp = gpar(col = "grey90", lwd=1, lty=2))
  grid.segments(c(1:(nmotif-1)/nmotif, rep(0, nmotif-1)), c((nmotif-1):1/nmotif, 1:(nmotif-1)/nmotif), c(1:(nmotif-1)/nmotif, (nmotif-1):1/nmotif), c(rep(1, nmotif-1), 1:(nmotif-1)/nmotif), gp = gpar(col = "grey90", lwd=2))

  vp <- viewport(layout=grid.layout(nmotif, nmotif)) #parse case
  pushViewport(vp)
  par(mar=c(2,1,1,1))
  for (line in 1:nmotif)
  {
    #plot distance
    if (line>1)
    {
      for (distance.col in 1:(line-1))
      {	
        #names
        if (bysim)
        {
          motifname1 <- similar(pos)[distance.col]
          motifname2 <- similar(pos)[line]					
        } 
        else 
        {
          motifname1 <- pos[[distance.col]]@motifName
          motifname2 <- pos[[line]]@motifName
        }

        vpdistance <- viewport(layout.pos.col=distance.col, layout.pos.row=line) #seqLogo
        pushViewport(vpdistance)	
        pushViewport(plotViewport(c(7/nmotif, 0.5,7/nmotif,0.5), name="plot"))
        if (is.null(xlim))
        {xlim <- c(-meanLength/2, meanLength/2)}
        xmin=xlim[1]
        xlength=xlim[2]-xlim[1]
        commonSequences <- unique(pos[[distance.col]]@positionVector$seq [pos[[distance.col]]@positionVector$seq %in% pos[[line]]@positionVector$seq]) #select common sequences 				
        if (!strand)
        {						
          dist <- motivDistance(pos[[distance.col]]@positionVector, pos[[line]]@positionVector, commonSequences, method, motifname1, motifname2, c("+","-"), meanLength)
          distance <- as.numeric(dist[[1]])
          distance.names <- dist[[2]]
          distance <- as.numeric(distance[distance>xlim[1] & distance<xlim[2]])

          if(length(distance[distance>xlim[1] & distance<xlim[2]])>1)
          {	#computes density
            dens <- density(as.numeric(distance[distance>xlim[1] & distance<xlim[2]]), bw=bw)
            density.x <- dens$x[dens$x>xlim[1] & dens$x<xlim[2]]
            density.y <- dens$y[dens$x>xlim[1] & dens$x<xlim[2]]
            x <- density.x - xmin
            y <- density.y
            ymax <- max(y)
            panel.lines(x=x/xlength, y=y/ymax, col="black", lwd=2)
            popViewport() #end plotViewport "plot"
            pushViewport(plotViewport(c(0, 0.5,0,0.5), name="axis"))
            panel.segments(0,0.13,1,0.13, col="black", lwd=1)
          }
          else
          {
            panel.text(x=0.5, y=0.5, labels="Not enough points", cex=2.5/nmotif)
            xlength=0
            xmin=0
          }
        }
        else if (strand)
        {	
          distp <- motivDistance(pos[[distance.col]]@positionVector, pos[[line]]@positionVector, commonSequences, method, motifname1, motifname2, "+", meanLength)
          distm <- motivDistance(pos[[distance.col]]@positionVector, pos[[line]]@positionVector, commonSequences, method, motifname1, motifname2, "-", meanLength)
          distance.positive <- distp[[1]]
          distance.negative <-	distm[[1]]				
          distance.names <- distp[[2]]
          if (length(distance.positive[distance.positive>xlim[1] & distance.positive<xlim[2]])>1 && length(distance.negative[distance.negative>xlim[1] & distance.negative<xlim[2]])>1) #both strands
          {
            density.positive <- density(distance.positive[distance.positive>xlim[1] & distance.positive<xlim[2]], bw=bw)
            density.positive.x <- density.positive$x[density.positive$x>xlim[1] & density.positive$x<xlim[2]]
			density.positive.y <- density.positive$y[density.positive$x>xlim[1] & density.positive$x<xlim[2]]
            density.negative <- density(distance.negative[distance.negative>xlim[1] & distance.negative<xlim[2]], bw=bw)
            density.negative.x <- density.negative$x[density.negative$x>xlim[1] & density.negative$x<xlim[2]]
            density.negative.y <- density.negative$y[density.negative$x>xlim[1] & density.negative$x<xlim[2]]
            density.positive.x.corrected <- density.positive.x - xmin
            density.positive.y.corrected <- density.positive.y
            density.negative.x.corrected <- density.negative.x- xmin
            density.negative.y.corrected <- density.negative.y		
            ymax <- max(density.positive.y.corrected, density.negative.y.corrected)	
			panel.lines(x=density.positive.x.corrected/xlength, y=density.positive.y.corrected/ymax, col="blue", lwd=2)
            panel.lines(x=density.negative.x.corrected/xlength, y=density.negative.y.corrected/ymax, col="red", lwd=2)						
          }
          else if (length(distance.positive[distance.positive>xlim[1] & distance.positive<xlim[2]])>1)
          { #forward strand only
            density.positive <- density(distance.positive, bw=bw)
            density.positive.x <- density.positive$x[density.positive$x>xlim[1] & density.positive$x<xlim[2]]
            density.positive.y <- density.positive$y[density.positive$x>xlim[1] & density.positive$x<xlim[2]]
            density.positive.x.corrected <- density.positive.x - xmin
            density.positive.y.corrected <- density.positive.y
            ymax <- max(density.positive.y.corrected)
            panel.lines(x=density.positive.x.corrected/xlength, y=density.positive.y.corrected/ymax, col="blue", lwd=2)							
          }
          else if (length(distance.negative[distance.negative>xlim[1] & distance.negative<xlim[2]])>1)
          { #reverse strand only
            density.negative <- density(distance.negative, bw=bw)
            density.negative.x <- density.negative$x[density.negative$x>xlim[1] & density.negative$x<xlim[2]]
            density.negative.y <- density.negative$y[density.negative$x>xlim[1] & density.negative$x<xlim[2]]
            density.positive.x.corrected <- density.negative.x - xmin
            density.positive.y.corrected <- density.negative.y
            ymax <- max(density.negative.y.corrected)
            panel.lines(x=density.negative.x.corrected/xlength, y=density.negative.y.corrected/ymax, col="red", lwd=2)			
          } 
          else
          {
            panel.text(x=0.5, y=0.5, labels="Not enough points", cex=2.5/nmotif)
            xlength=0
            xmin=0
          }
        }
        popViewport() #end plotViewport "plot"
        pushViewport(plotViewport(c(0, 0.5,0,0.5), name="axis"))
        panel.segments(0,0.13,1,0.13, col="black", lwd=1)

        x.zero=(0-xmin)/xlength
        x.firstQuarter=x.zero*1/2
        x.thirdQuarter=(1-x.zero)/2+x.zero
        panel.segments(x0=c(x.zero, x.firstQuarter, x.thirdQuarter), y0=c(0.1, 0.1, 0.1), x1=c(x.zero, x.firstQuarter, x.thirdQuarter), y1=c(0.155, 0.155, 0.155), col="black", lwd=1)	
        panel.text(x=x.zero, y=0.045, labels="0", cex=2.5/nmotif)
        panel.text(x=x.firstQuarter, y=0.045, labels=round(x.firstQuarter*xlength+xmin), cex=2.5/nmotif)
        panel.text(x=x.thirdQuarter, y=0.045, labels=round(x.thirdQuarter*xlength+xmin), cex=2.5/nmotif)
        popViewport() #end plotViewport "axis"
        grid.text(x=unit(0.5,"npc"), y=unit(0.92,"npc"), paste("d(", distance.names[1],"-", distance.names[2],")"))
        popViewport() #end vpdistance
      }
    }

    #print motif name
    vpname <- viewport(layout.pos.col=line, layout.pos.row=line) #seqLogo
    pushViewport(vpname)	
    pushViewport(	viewport(x=0.5, y=0.5, width=0.6, height=0.25, angle=45, name="rect"))
    grid.rect(x=unit(0.5,"npc"), y=unit(0.5,"npc"), height=unit(1,"npc"), width=unit(1,"npc"), gp=gpar(col="grey95", fill="grey95"))				

    if (bysim)
    {
      grid.text(x=unit(0.5,"npc"), y=unit(0.5,"npc"), pos[[line]]@similarMotif, gp=gpar(cex=6/nmotif, srt=45, col="red", font=2))	
    } 
    else 
    {
      grid.text(x=unit(0.5,"npc"), y=unit(0.5,"npc"), pos[[line]]@motifName, gp=gpar(cex=6/nmotif, srt=45, col="red", font=2))	
    }
    popViewport()	 #end "rect"
    popViewport()	#end vpname

    #plot common sequences
    if (line<nmotif)
    {
      for(venn.col in (line+1):nmotif)
      {
        vpseqcom <- viewport(layout.pos.col=venn.col, layout.pos.row=line) #seqLogo
        pushViewport(vpseqcom)	
        numberOfCommonSequences=length(unique(pos[[line]]@positionVector$seq [ pos[[line]]@positionVector$seq %in% pos [[venn.col]]@positionVector$seq] ) )#compte nombre sequences communes
        numberOfSequences1=length(unique(pos[[line]]@positionVector$seq))
        numberOfSequences2=length(unique(pos[[venn.col]]@positionVector$seq))
        numberOfSequences1.single=numberOfSequences1-numberOfCommonSequences
        numberOfSequences2.single=numberOfSequences2-numberOfCommonSequences

        grid.circle(x=unit(0.4,"npc"), y=unit(0.5,"npc"), r=unit(0.3,"npc"), gp=gpar(col="blue", lwd=2))
        grid.circle(x=unit(0.6,"npc"), y=unit(0.5,"npc"), r=unit(0.3,"npc"), gp=gpar(col="green", lwd=2))
        if (bysim)
        {
          grid.text(x=unit(0.1,"npc"), y=unit(0.9,"npc"), pos[[line]]@similarMotif, just="left", gp=gpar(col="blue"))
          grid.text(x=unit(0.9,"npc"), y=unit(0.9,"npc"), pos[[venn.col]]@similarMotif, just="right", gp=gpar(col="green"))
        }
        else
        {
          grid.text(x=unit(0.1,"npc"), y=unit(0.9,"npc"), pos[[line]]@motifName, just="left", gp=gpar(col="blue"))
          grid.text(x=unit(0.9,"npc"), y=unit(0.9,"npc"), pos[[venn.col]]@motifName, just="right", gp=gpar(col="green"))
        }
        grid.text(x=unit(0.5,"npc"), y=unit(0.5,"npc"), numberOfCommonSequences, just="center", gp=gpar(col="black"))
        grid.text(x=unit(c(0.35,0.35),"npc"), y=unit(c(0.5,0.35),"npc"), c(numberOfSequences1.single, paste(round(100*numberOfSequences1.single/numberOfSequences1,1),"%", sep="")), just="center", gp=gpar(col="blue", cex=c(1,0.7)))
        grid.text(x=unit(c(0.65,0.65),"npc"), y=unit(c(0.5,0.35),"npc"), c(numberOfSequences2.single, paste(round(100*numberOfSequences2.single/numberOfSequences2,1),"%", sep="")), just="center", gp=gpar(col="green", cex=c(1,0.7)))
        popViewport() #end vpseqcom	
      }
    }
  }
  popViewport() #end vp
  popViewport() #end plotViewport "grid"
}
