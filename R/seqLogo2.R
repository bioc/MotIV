
# SeqLogo
# Author  	Oliver Bembom
# Maintainer 	Oliver Bembom

#modified by Eloi Mercier

######################
##
## plot sequence logo
##
######################

letterA <- function(x.pos, y.pos, ht, wt, id=NULL)
{
  x <- c(0,4,6,10,8,6.8,3.2,2,0,3.6,5,6.4,3.6)
  y <- c(0,10,10,0,0,3,3,0,0,4,7.5,4,4)
  x <- 0.1*x
  y <- 0.1*y

  x <- x.pos + wt*x
  y <- y.pos + ht*y

  if (is.null(id))
  {
    id <- c(rep(1,9), rep(2,4))
  }
  else
  {
    id <- c(rep(id,9), rep(id+1,4))
  }

  fill <- c("green","white")

  list(x=x, y=y, id=id, fill=fill)
}

## T
letterT <- function(x.pos, y.pos, ht, wt, id=NULL)
{
  x <- c(0,10,10,6,6,4,4,0)
  y <- c(10,10,9,9,0,0,9,9)
  x <- 0.1*x
  y <- 0.1*y

  x <- x.pos + wt*x
  y <- y.pos + ht*y

  if (is.null(id)){
    id <- rep(1,8)
  }
  else
  {
    id <- rep(id,8)
  }

  fill <- "red"

  list(x=x, y=y, id=id, fill=fill)
}

## C
letterC <- function(x.pos, y.pos, ht, wt, id=NULL)
{
  angle1 <- seq(0.3+pi/2, pi, length=100)
  angle2 <- seq(pi,1.5*pi, length=100)
  x.l1 <- 0.5 + 0.5*sin(angle1)
  y.l1 <- 0.5 + 0.5*cos(angle1)
  x.l2 <- 0.5 + 0.5*sin(angle2)
  y.l2 <- 0.5 + 0.5*cos(angle2)

  x.l <- c(x.l1, x.l2)
  y.l <- c(y.l1, y.l2)

  x <- c(x.l, rev(x.l))
  y <- c(y.l,1-rev(y.l))

  x.i1 <- 0.5 +0.35*sin(angle1)
  y.i1 <- 0.5 +0.35*cos(angle1)
  x.i1 <- x.i1[y.i1<=max(y.l1)]
  y.i1 <- y.i1[y.i1<=max(y.l1)]
  y.i1[1] <- max(y.l1)

  x.i2 <- 0.5 +0.35*sin(angle2)
  y.i2 <- 0.5 +0.35*cos(angle2)

  x.i <- c(x.i1, x.i2)
  y.i <- c(y.i1, y.i2)

  x1 <- c(x.i, rev(x.i))
  y1 <- c(y.i,1-rev(y.i))

  x <- c(x, rev(x1))
  y <- c(y, rev(y1))

  x <- x.pos + wt*x
  y <- y.pos + ht*y

  if (is.null(id))
  {
    id <- rep(1, length(x))
  }
  else
  {
    id <- rep(id, length(x))
  }

  fill <- "blue"

  list(x=x, y=y, id=id, fill=fill)
}

## G
letterG <- function(x.pos, y.pos, ht, wt, id=NULL)
{
  angle1 <- seq(0.3+pi/2, pi, length=100)
  angle2 <- seq(pi,1.5*pi, length=100)
  x.l1 <- 0.5 + 0.5*sin(angle1)
  y.l1 <- 0.5 + 0.5*cos(angle1)
  x.l2 <- 0.5 + 0.5*sin(angle2)
  y.l2 <- 0.5 + 0.5*cos(angle2)
  
  x.l <- c(x.l1, x.l2)
  y.l <- c(y.l1, y.l2)
  
  x <- c(x.l, rev(x.l))
  y <- c(y.l,1-rev(y.l))
  
  x.i1 <- 0.5 +0.35*sin(angle1)
  y.i1 <- 0.5 +0.35*cos(angle1)
  x.i1 <- x.i1[y.i1<=max(y.l1)]
  y.i1 <- y.i1[y.i1<=max(y.l1)]
  y.i1[1] <- max(y.l1)
  
  x.i2 <- 0.5 +0.35*sin(angle2)
  y.i2 <- 0.5 +0.35*cos(angle2)
  
  x.i <- c(x.i1, x.i2)
  y.i <- c(y.i1, y.i2)
  
  x1 <- c(x.i, rev(x.i))
  y1 <- c(y.i,1-rev(y.i))
  
  x <- c(x, rev(x1))
  y <- c(y, rev(y1))
  
  h1 <- max(y.l1)
  r1 <- max(x.l1)
  
  h1 <- 0.4
  x.add <- c(r1,0.5,0.5, r1-0.2, r1-0.2, r1, r1)
  y.add <- c(h1, h1, h1-0.1, h1-0.1,0,0, h1)

  

  if (is.null(id))
  {
    id <- c(rep(1, length(x)), rep(2, length(x.add)))
  }
  else
  {
    id <- c(rep(id, length(x)), rep(id+1, length(x.add)))
  }

  x <- c(rev(x), x.add)
  y <- c(rev(y), y.add)
    
  x <- x.pos + wt*x
  y <- y.pos + ht*y
  
  
  fill <- c("orange","orange")
    
  list(x=x, y=y, id=id, fill=fill)
}


addLetter <- function(letters, which, x.pos, y.pos, ht, wt)
{
  if (which == "A")
  {
    letter <- letterA(x.pos, y.pos, ht, wt)
  }
  else if (which == "C")
  {
    letter <- letterC(x.pos, y.pos, ht, wt)    
  }
  else if (which == "G")
  {
    letter <- letterG(x.pos, y.pos, ht, wt)    
  }
  else if (which == "T")
  {
    letter <- letterT(x.pos, y.pos, ht, wt)    
  }
  else
  {
    stop("which must be one of A, C, G, T")
  }

  letters$x <- c(letters$x, letter$x)
  letters$y <- c(letters$y, letter$y)

  lastID <- ifelse(is.null(letters$id),0, max(letters$id))
  letters$id <- c(letters$id, lastID+letter$id)
  letters$fill <- c(letters$fill, letter$fill)
  letters
}

###############
makePWM <- function(pwm, alphabet="DNA")
{

  if (is.data.frame(pwm)) pwm <- as.matrix(pwm)
  if (!is.matrix(pwm)) stop("pwm must be a matrix or a dataframe")

  if (!alphabet %in% c("DNA","AA"))
  {
    stop("alphabet must be either DNA or AA")
  }
  if (alphabet == "DNA" & nrow(pwm) != 4)
  {
    stop("PWM for DNA motifs must have 4 rows")
  }
  else if (alphabet == "AA" && nrow(pwm) != 21)
  {
    stop("PWM for amino acid motifs must have 21 rows")
  }
  if (any(abs(1 - apply(pwm,2, sum)) > 0.01))
  {
    print(apply(pwm,2, sum))
    warning("Columns of PWM must add up to 1.0")
  }

  width <- ncol(pwm)
  colnames(pwm) <- 1:width
  rownames(pwm) <- c("A","C","G","T")

  cons <- pwm2cons(pwm)
  ic <- pwm2ic(pwm)

  new("pwm2", pwm=pwm, consensus=cons, ic=ic, width=width, alphabet=alphabet)
}



## get information content profile from PWM
pwm2ic <- function(pwm)
{
    npos <- ncol(pwm)
    ic <- numeric(length=npos)
    for (i in 1:npos) {
        ic[i]<-2 + sum(sapply(pwm[, i], function(x) { if (x > 0){x*log2(x);}else{0;}}))
    }
    ic
}

## get consensus sequence from PWM
pwm2cons <- function(pwm)
{
    if (class(pwm)!="matrix") {warning("pwm argument must be of class matrix")}
    letters <- c("A", "C", "G", "T")
    paste(apply(pwm, 2, function(x){letters[rev(order(x))[1]]}), collapse="")
}


## plot a sequence logo

seqLogo2 <- function(pwm, ic.scale=TRUE, xaxis=TRUE, yaxis=TRUE, xfontsize=15, yfontsize=15, vmargins=c(0,0), hmargins=c(0,0), size=1, trim=0)
{
  if (class(pwm) == "pwm2")
  {
    pwm <- pwm@pwm    
  }
  else if (class(pwm) == "data.frame")
  {
    pwm <- as.matrix(pwm)
  }
  else if (class(pwm) != "matrix")
  {
    stop("pwm must be of class matrix or data.frame")
  }

  if (any(abs(1 - apply(pwm,2, sum)) > 0.01))
  stop("Columns of PWM must add up to 1.0")


  chars <- c("A","C","G","T")
  letters <- list(x=NULL, y=NULL, id=NULL, fill=NULL)
  npos <- ncol(pwm)


  if (ic.scale)
  {
    ylim <-2*size
    ylab <- "Information content"
    facs <- pwm2ic(pwm)
  }
  else
  {
    ylim <- 1*size
    ylab <- "Probability"
    facs <- rep(1, npos)
  }

  wt <- 1  
  x.pos <- 0  
  for (j in 1:npos)
  {

    column <- pwm[, j]
    hts <- column*facs[j]*size
    letterOrder <- order(hts)

    y.pos <- 0    
    for (i in 1:4)
    {
      letter <- chars[letterOrder[i]]
      ht <- hts[letterOrder[i]]
      if (ht>trim) 
	 {
	 letters <- addLetter(letters, letter, x.pos, y.pos, ht, wt)
     y.pos <- y.pos + ht + 0.01
	 }
    }
    x.pos <- x.pos + wt
  }

  #grid.newpage()

  pushViewport(plotViewport(c(vmargins[1] + xfontsize/(3.5*size), hmargins[1]+ xfontsize/(3.5*size) , vmargins[2], hmargins[2] ))) #margins
  pushViewport(dataViewport(0:ncol(pwm)/size,0:ylim/size, name="vp1"))
  grid.polygon(x=unit(letters$x,"native"), y=unit(letters$y,"native"),
  id=letters$id,
  gp=gpar(fill=letters$fill, col="transparent"))
  if (xaxis)
  {
    grid.xaxis(at=seq(0.5, ncol(pwm)-0.5), label=1:ncol(pwm), gp=gpar(fontsize=xfontsize*size))
    grid.text("Position", x=unit(0.5*size,"snpc"), y=unit(-0.1*size,"npc"), gp=gpar(fontsize=xfontsize*size))
  }
  if (yaxis)
  {
    grid.yaxis(at=seq(from=0, to=ylim, by=0.5), label=seq(from=0, to=ylim, by=0.5), gp=gpar(fontsize=yfontsize*size))
    grid.text(ylab, x=unit(-1,"lines"), rot=90, gp=gpar(fontsize=yfontsize*size))
  }
  popViewport()
  popViewport()
  par(ask=FALSE)
}
