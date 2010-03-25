.moduleRanges <- function(lx,ly, range)
{
  data <- NULL
  for(i in 1:length(lx))
  {
    for (j in 1:length(ly))
    {
      if(abs(lx[i]-ly[j])<=range)
      {
        data$x=c(data$x,i)
        data$y=c(data$y,j)
      }
    }
  }
  return(data)
}

################


calculateModules <- function(x, y, range=250, merge=TRUE) # x and y should be RangedData objects
{
  iranges <- NULL
  data <- RangedData()
  chr.common <- which(unique(x[["chr"]]) %in% unique(y[["chr"]]))
  chr.common.name=unique(x[["chr"]])[chr.common]

  for (i in chr.common)
  {
    chr.crt=unique(x[["chr"]])[i]
    chr1.x=which(x[["chr"]]==chr.crt)
    chr1.y=which(y[["chr"]]==chr.crt)
    xranges=ranges(x)[[1]][chr1.x]
    yranges=ranges(y)[[1]][chr1.y]

    start <- NULL
    end <- NULL
    motif1 <- NULL
    motif2 <- NULL

    modules <- .moduleRanges(start(xranges), start(yranges), range) #selects close motifs
    if (merge)
    { 
      for (j in seq(length(modules$x)))
      {
        if (min(xranges[modules$x[j]]) < min(yranges[modules$y[j]]))
        {
          start.min = min(xranges[modules$x[j]])
          end.max = max(yranges[modules$y[j]])
          start.motif = x[["name"]][1]
          end.motif = y[["name"]][1]
        } 
        else
        {
          start.min = min(yranges[modules$y[j]])
          end.max = max(xranges[modules$x[j]])
          start.motif = y[["name"]][1]
          end.motif = x[["name"]][1]
        }

        start <- c(start, start.min)
        end <- c(end, end.max)
        motif1 <- c(motif1, start.motif)
        motif2 <- c(motif2, end.motif)	
      }

      chr <- rep(chr.crt,length(start))
      chr.strand <- rep("+",length(start))
      iranges <- IRanges(start=start,end=end)
      data <- rbind( data, RangedData( iranges, chr,strand=chr.strand, motif1=motif1, motif2=motif2))
    }
    else
    {
      data <- rbind(data, x[modules$x,], y[modules$y,])
    }
  }
  return(data)
}