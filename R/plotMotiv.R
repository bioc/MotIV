

plotMotif <- function( object, current, column,  top, layout, bysim, group, trim, cex)
{
	size=max(layout[1], layout[2])		
	
	if (group || bysim)
	{
		plotname <- similarity(object)
	}         else         {
		plotname <-switch (class(object), "motiv"=names(object), "list"=sapply(object, function(x){x@motifName}))  
	}

	pushViewport(plotViewport(c(0, 0,2/layout[1],0)))
	vpin <- viewport(layout=grid.layout(top+1,4)) #parse case
	pushViewport(vpin)
	
	if (top > 1 && class(object)[[1]]=="motiv")
	{
		grid.segments( x0=rep(0, top), y0=c(0:(top-1)/(top+1)), x1=rep(1, top), y1=c(0:(top-1)/(top+1)), gp = gpar(col = "grey90", lwd=1))
	}
	
#plot motif and reverse
	mat <- switch (class(object), "motiv"=object@input[[current]], "list"=object[[current]]@pwm[[1]])
	matrix <- apply(mat,2, function(mat, column){mat[column]/sum(mat[column])})
	revmat <- matrix[4:1, dim(matrix)[2]:1]
	vpmotif <- viewport(layout.pos.col=1:2, layout.pos.row=1) #seqLogo
	pushViewport(vpmotif)	
	grid.rect(x=unit(1,"npc"), y=unit(0.45,"npc"), height=unit(0.95,"npc"), width=unit(1.9,"npc"), gp=gpar(col="grey95", fill="grey95"))	
	seqLogo2(matrix, yaxis=F, xaxis=F, size=1, yfontsize=0, xfontsize=0, vmargins=c(0,0.4), hmargins=c(0.1,0.1), trim=trim)
	popViewport() #end vpmotif
	vpmotifrev <- viewport(layout.pos.col=3:4, layout.pos.row=1) #seqLogo
	pushViewport(vpmotifrev)	
	seqLogo2(revmat, yaxis=F, xaxis=F, size=1, yfontsize=0, xfontsize=0, vmargins=c(0, 0.4), hmargins=c(0.1, 0.1))
	
	do.call("grid.text", list(c(plotname[current],"forward","RC"), x=unit(c(0,-0.5,0.5),"npc"), y=unit(c(1.05,1,1),"npc"), gp=gpar(cex=c(max(3/size,0.7), rep(max(2/size,0.6),2))*cex, font=c(2,3,3))))	
}

#####MOTIV#####

plotMotiv <- function (motiv, ncol, nrow, top, bysim, rev, main, sub, trim, cex)
{
	if(ncol==0 && nrow==0)
	{
		ncol=ceiling(sqrt(length(motiv@bestMatch)))
		nrow=ncol
	} 	else if (ncol==0 && nrow!=0)	{
		ncol=ceiling(length(motiv@bestMatch)/nrow)
	}	else if ((ncol!=0 && nrow==0))	{
		nrow=ceiling(length(motiv@bestMatch)/ncol)
	}
	
	layout <- c(nrow, ncol)
	if (top > length(motiv@bestMatch[[1]]@aligns) || top <1)
	{top <- length(motiv@bestMatch[[1]]@aligns )  }
	inputPWM <- motiv@input	
	size=max(layout[1], layout[2])		
	
#plot grid and title
	grid.newpage()
	grid.text(sub, gp=gpar(col="black", font=3, cex=1.3*cex), y=unit(0.03,"npc"), just="top")
	
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
				
				plotMotif( motiv, p, i,  top, layout, bysim, FALSE, trim=trim, cex=cex)
				popViewport() #end vpmotifrev
				
#plot best matches
				for (k in 1:top)
				{
					strand <- "+"	
					vptop3 <- viewport(layout.pos.col=1:2, layout.pos.row=k+1)
					pushViewport(vptop3)		
					mat3<- motiv@bestMatch[[p]]@aligns[[k]]@TF@pwm
					if (motiv@bestMatch[[p]]@aligns[[k]]@strand=="-") 
					{
						if (rev)
						{
							mat3=mat3[4:1, dim(mat3)[2]:1]
						}
						else
						{
							strand <- "-" 
						}
					}
					mat3<- makePWM(mat3)
					seqLogo2(mat3, yaxis=F, xaxis=F, size=1, yfontsize=0, xfontsize=0, vmargins=c(0.1, 0.1), hmargins=c(0.2,0.2), trim=trim) #plot Logo
					popViewport() #end vptop3
					vpeval <- viewport(layout.pos.col=3:4, layout.pos.row=k+1) 
					pushViewport(vpeval)		
					grid.text(paste(motiv@bestMatch[[p]]@aligns[[k]]@TF@name), gp=gpar(col="black", cex=max(0.6,2/size)*cex, font=2), y=unit(0.75,"npc"))
					grid.text(paste(strand), gp=gpar(col="black", cex=max(0.6,2.2/size)*cex), y=unit(0.5,"npc"))
					grid.text(paste( format(motiv@bestMatch[[p]]@aligns[[k]]@evalue, scientific=T, digits=5)), gp=gpar(col="black", cex=max(0.6,2/size)*cex), y=unit(0.2,"npc"))
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
		grid.text(main, gp=gpar(col="black", font=2, cex=max(0.8,3/size)), y=unit(0.999,"npc"), just="top")

}