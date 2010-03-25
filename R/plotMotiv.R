
#####MOTIV#####

plotMotiv <- function (motiv, ncol, nrow, top, bysim, rev, main, sub)
{
  layout <- c(nrow, ncol)
  inputPWM <- motiv@input	
  size=max(layout[1], layout[2])		
  similarMotifs=1:length(similar(motiv))
  names(similarMotifs)=similar(motiv)

  #plot grid and title
  grid.newpage()
  grid.text(main, gp=gpar(col="black", font=2, cex=1.5), y=unit(0.98,"npc"))
  grid.text(sub, gp=gpar(col="black", font=3, cex=1.3), y=unit(0.03,"npc"), just="top")
  grid.text("RC : Reverse Complement", gp=gpar(col="black", font=3, cex=0.7), y=unit(0.015,"npc"), x=unit(0.99,"npc"), just="right")
  if (!bysim)
  {
    grid.text("* : similar motifs", gp=gpar(col="black", font=3, cex=0.7), y=unit(0.015,"npc"), x=unit(0.01,"npc"), just="left")
  }

  pushViewport(plotViewport(c(7/layout[1], 0,7/layout[1],0)))
  vp  <- viewport(layout=grid.layout(layout[1], layout[2]))	#1rst lvl
  pushViewport(vp)	
  grid.segments(c(1:(layout[2]-1)/layout[2], rep(0, layout[1]-1)), c(rep(0, layout[2]-1), 1:(layout[1]-1)/layout[1]), c(1:(layout[2]-1)/layout[2], rep(1, layout[1]-1)), c(rep(1, layout[2]-1), 1:(layout[1]-1)/layout[1]), gp = gpar(col = "grey90", lwd=2))
  #plot Logo	
  p=1
  for (j in 1:layout[1])
  {
    for (i in 1:layout[2])
    {
      if (p<=length(inputPWM))
      {
        vpcase <- viewport(layout.pos.col=i, layout.pos.row=j)
        pushViewport(vpcase)	
        #names
        if (bysim)
        {
          plotname <- motiv@bestMatch[[p]]@similar
        } 
        else 
        {
          n=similarMotifs[motiv@bestMatch[[p]]@similar]
          plotname <- paste(c(motiv@bestMatch[[p]]@name , rep("*", n)), sep="", collapse="")
        }

        grid.text(c(plotname,"forward","RC"), x=unit(c(0.5,0.25,0.75),"npc"), y=unit(0.96,"npc"), gp=gpar(cex=c(3/size, rep(max(2/size,0.6),2)), font=c(2,3,3)))	
        pushViewport(plotViewport(c(0, 0,2/layout[1],0)))
        vpin <- viewport(layout=grid.layout(top+1,4)) #parse case
        pushViewport(vpin)

        if (top > 1)
        {
          grid.segments( x0=rep(0, top), y0=c(0:(top-1)/(top+1)), x1=rep(1, top), y1=c(0:(top-1)/(top+1)), gp = gpar(col = "grey90", lwd=1))
        }

        #plot motif and reverse
        mat <- inputPWM[[p]]
        mat <- apply(mat,2, function(mat, i){mat[i]/sum(mat[i])})
        revmat <- mat[4:1, dim(mat)[2]:1]
        vpmotif <- viewport(layout.pos.col=1:2, layout.pos.row=1) #seqLogo
        pushViewport(vpmotif)	
        grid.rect(x=unit(1,"npc"), y=unit(0.45,"npc"), height=unit(0.95,"npc"), width=unit(1.9,"npc"), gp=gpar(col="grey95", fill="grey95"))	
        seqLogo(mat, yaxis=F, xaxis=F, size=1, yfontsize=0, xfontsize=0, vmargins=0.1, hmargins=0.1)
        popViewport() #end vpmotif
        vpmotifrev <- viewport(layout.pos.col=3:4, layout.pos.row=1) #seqLogo
        pushViewport(vpmotifrev)	
        seqLogo(revmat, yaxis=F, xaxis=F, size=1, yfontsize=0, xfontsize=0, vmargins=0.1, hmargins=0.1)
        grid.text("Best Matches", x=unit(0.5,"npc"), y=unit(-0.15,"npc"), gp=gpar(cex=max(0.6,2.5/size), font=2))	
        popViewport() #end vpmotifrev

        #plot best matches
        for (top3 in 1:top)
        {
          rc <- NULL
          vptop3 <- viewport(layout.pos.col=1:2, layout.pos.row=top3+1)
          pushViewport(vptop3)		
          mat3<- motiv@bestMatch[[p]]@aligns[[top3]]@TF@pwm
          if (motiv@bestMatch[[p]]@aligns[[top3]]@strand=="-") 
          {
            if (rev)
            {
              mat3=mat3[4:1, dim(mat3)[2]:1]
            }
            else
            {
              rc="(RC)"
            }
          }
          mat3<- makePWM(mat3)
          seqLogo(mat3, yaxis=F, xaxis=F, size=1, yfontsize=0, xfontsize=0, vmargins=0.1, hmargins=0.2) #plot Logo
          popViewport() #end vptop3
          vpeval <- viewport(layout.pos.col=3:4, layout.pos.row=top3+1) 
          pushViewport(vpeval)		
          grid.text(paste(motiv@bestMatch[[p]]@aligns[[top3]]@TF@name, rc), gp=gpar(col="black", cex=max(0.6,2/size)), y=unit(0.55,"npc"))
          grid.text(paste("Evalue : ", format(motiv@bestMatch[[p]]@aligns[[top3]]@evalue, scientific=T, digits=5)), gp=gpar(col="black", cex=max(0.6,2/size)), y=unit(0.2,"npc"))
          popViewport()	#end vpeval	
        }
        popViewport() #end vpin
        popViewport() #end vpcase
        popViewport() #end plotViewport "case"
        p=p+1
      }
    }
  }
  popViewport() #end vp
  popViewport() #end plotViewport "grid"
}