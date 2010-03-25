
#####DISTRIBUTION#####
plotDistribution <- function (pos, group, main, sort, ncol, nrow, strand, bysim, meanLength)
{
  positionVector <- list()
  for (i in 1:length(pos))
  {
    positionVector[[i]]=(pos[[i]]@positionVector$end+ pos[[i]]@positionVector$start)/2
  }
  if (group) 
  {
    names(positionVector)=similar(pos)
  } 
  else 
  {
    names(positionVector)=sapply(pos, function(x){x@motifName})
  }
  var <- sapply(positionVector, function(x){var(as.integer(x))})
  var.reference=which(var==min(var))
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
  similarMotifs=1:length(unique(similar(pos)))
  names(similarMotifs)=unique(similar(pos))

  grid.newpage()
  grid.text(main, gp=gpar(col="black", font=2, cex=1.5), y=unit(0.98,"npc"))
  grid.text("position", gp=gpar(col="black", font=3, cex=1), y=unit(0.01,"npc"))
  grid.text("RC : Reverse Complement", gp=gpar(col="black", font=3, cex=0.7), y=unit(0.015,"npc"), x=unit(0.99,"npc"), just="right")
  if (!bysim)
  {
    grid.text("* : similar motifs", gp=gpar(col="black", font=3, cex=0.7), y=unit(0.015,"npc"), x=unit(0.01,"npc"), just="left")
  }

  if (strand)
  {
    grid.text("Forward", gp=gpar( font=3, cex=0.7,  col="blue"), y=unit(0.02,"npc"), x=unit(0.18,"npc"), just="top")
    grid.text("Backward", gp=gpar( font=3, cex=0.7,  col="red"), y=unit(0.02,"npc"), x=unit(0.27,"npc"), just="top")
  }
  pushViewport(plotViewport(c(7/layout[1], 0,7/layout[1],0), name="grid"))
  vp  <- viewport(layout=grid.layout(layout[1], layout[2]))	#1rst lvl
  pushViewport(vp)
  grid.segments(c(1:(layout[2]-1)/layout[2], rep(0, layout[1]-1)), c(rep(0, layout[2]-1), 1:(layout[1]-1)/layout[1]), c(1:(layout[2]-1)/layout[2], rep(1, layout[1]-1)), c(rep(1, layout[2]-1), 1:(layout[1]-1)/layout[1]), gp = gpar(col = "grey90", lwd=2))

  #plot Logo
  countMotifs=1
  for (line in 1:layout[1])
  {
    for (colum in 1:layout[2])
    {
      if (countMotifs<=length(pos))
      {
        curr=which(names(var.sort)[countMotifs]==names(var))
        if (group || bysim)
        {
          plotname <- pos[[curr]]@similarMotif
        }
        else
        {
          n=similarMotifs[pos[[curr]]@similarMotif]
          plotname <- paste(c(pos[[curr]]@motifName , rep("*", n)), sep="", collapse="")
        }
        vpcase <- viewport(layout.pos.col=colum, layout.pos.row=line)
        pushViewport(vpcase)
        grid.text(c(plotname,"forward","RC"), x=unit(c(0.5,0.25,0.75),"npc"), y=unit(0.96,"npc"), gp=gpar(cex=c(3/size, rep(max(2/size,0.6),2)), font=c(2,3,3)))
        pushViewport(plotViewport(c(0, 0,2/layout[1],0), name="case"))
        vpin <- viewport(layout=grid.layout(4,4)) #parse case
        pushViewport(vpin)
        #plot motif and reverse
        mat <- pos[[curr]]@pwm[[1]]
        mat <- apply(mat,2, function(mat, colum){mat[colum]/sum(mat[colum])})
        revmat <- mat[4:1, dim(mat)[2]:1]
        vpmotif <- viewport(layout.pos.col=1:2, layout.pos.row=1) #seqLogo
        pushViewport(vpmotif)
        grid.rect(x=unit(1,"npc"), y=unit(0.45,"npc"), height=unit(0.95,"npc"), width=unit(1.9,"npc"), gp=gpar(col="grey95", fill="grey95"))
        seqLogo(mat, yaxis=F, xaxis=F, size=1, yfontsize=0, xfontsize=0, vmargins=0.1, hmargins=0.1)
        popViewport() #end vpmotif
        vpmotifrev <- viewport(layout.pos.col=3:4, layout.pos.row=1) #seqLogo
        pushViewport(vpmotifrev)
        seqLogo(revmat, yaxis=F, xaxis=F, size=1, yfontsize=0, xfontsize=0, vmargins=0.1, hmargins=0.1)
        grid.text("Distribution", x=unit(0,"npc"), y=unit(-0.15,"npc"), gp=gpar(cex=max(0.6,2.5/size), font=2))
        popViewport() #end vpmotifrev

        #plot distribution
        vpdistribution <- viewport(layout.pos.col=2:3, layout.pos.row=2:4) #seqLogo
        pushViewport(vpdistribution)		
        pushViewport(plotViewport(c(max(5/size,0.9), 0,3.5/size,0), name="plot"))

        if (!strand)
        {
          pos.start <- (pos[[curr]]@positionVector$end + pos[[curr]]@positionVector$start)/2
          dens <- density(pos.start )
          density.x <- dens$x - min(dens$x)
          density.y <- dens$y - min(dens$y)
          panel.lines(x=density.x/max(density.x), y=density.y/max(density.y), col="black", lwd=2)
          pos.start.center <- pos.start-mean(pos.start)
          kurtosis=mean(pos.start.center ^4)/mean(pos.start.center ^2)^2-3
          grid.text(paste("Kurtosis\n", round(kurtosis,3)), x=unit(-0.2,"npc"), y=unit(0.6,"npc"), gp=gpar(cex=max(0.6,2/size), font=2))
          grid.text(paste("sd\n", round(sd(pos.start),3)), x=unit(1.2,"npc"), y=unit(0.6,"npc"), gp=gpar(cex=max(0.6,2/size), font=2))
        } 
        else if (strand)
        {
          pos.start.positive <- (pos[[curr]]@positionVector$end[pos[[curr]]@positionVector$gadem.strand=="+"] + pos[[curr]]@positionVector$start[pos[[curr]]@positionVector$gadem.strand=="+"])/2 
          densp <- density(pos.start.positive)
          density.positive.x.corrected <- densp$x - min(densp$x)
          density.positive.y.corrected <- densp$y - min(densp$y)
          panel.lines(x=density.positive.x.corrected/max(density.positive.x.corrected), y=density.positive.y.corrected/max(density.positive.y.corrected), col="blue", lwd=2, lty=1)
          pos.start.positive.center <- pos.start.positive-mean(pos.start.positive)

          pos.start.negative <- (pos[[curr]]@positionVector$end[pos[[curr]]@positionVector$gadem.strand=="-"] + pos[[curr]]@positionVector$start[pos[[curr]]@positionVector$gadem.strand=="-"])/2
          densm <- density(pos.start.negative)
          density.negative.x.corrected <- densm$x - min(densm$x)
          density.negative.y.corrected <- densm$y - min(densm$y)
          panel.lines(x=density.negative.x.corrected/max(density.negative.x.corrected), y=density.negative.y.corrected/max(density.negative.y.corrected), col="red", lwd=2, lty=1)
          pos.start.negative.center <- pos.start.negative-mean(pos.start.negative)

          kurtosis.positive=mean(pos.start.positive.center ^4)/mean(pos.start.positive.center ^2)^2-3
          kurtosis.negative=mean(pos.start.negative.center^4)/mean(pos.start.negative.center^2)^2-3
          grid.text(paste("Kurtosis\n", round(kurtosis.positive,3)), x=unit(-0.2,"npc"), y=unit(0.6,"npc"), gp=gpar(cex=max(0.6,2/size), font=2, col="blue"))
          grid.text(paste(round(kurtosis.negative,3)), x=unit(-0.2,"npc"), y=unit(0.2,"npc"), gp=gpar(cex=max(0.6,2/size), font=2, col="red"))
          grid.text(paste("sd\n", round(sd(pos.start.positive),3)), x=unit(1.2,"npc"), y=unit(0.6,"npc"), gp=gpar(cex=max(0.6,2/size), font=2, col="blue"))
          grid.text(paste(round(sd(pos.start.negative),3)), x=unit(1.2,"npc"), y=unit(0.2,"npc"), gp=gpar(cex=max(0.6,2/size), font=2, col="red"))
        }
        popViewport() #end vpdistribution
        panel.segments(0,0.13,1,0.13, col="black", lwd=1)
        panel.segments(x0=c(0.25,0.5,0.75), y0=0.1+c(0,0,0), x1=c(0.25,0.5,0.75), y1=0.115+c(0.05,0.05,0.05), col="black", lwd=1)	
        panel.text(x=c(0.25,0.5,0.75), y=0.045, labels=round(c(-meanLength/4,0, meanLength/4),0), cex=2/size)
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
