
#####RANGED DATA#####
exportAsRangedData <- function (x, y, bysim=TRUE, correction=TRUE)
{
  if (class(x)!="motiv")
  {
    stop("x must be an object of class motiv.")
  }
  data <- RangedData()
  curr <- which(names(y)%in%names(x))
  vector.pos <- calculatePositionVector(x, y, group=F, correction=correction)
  p=1
  for (i in curr)
  {
    valid <- x@bestMatch[[p]]@valid 
    #names
    if (bysim)
    {
      mot <- rep(x@bestMatch[[p]]@similar, length(vector.pos[[p]]@positionVector$start))
    } 
    else 
    {
      mot <- rep(x@bestMatch[[p]]@name, length(vector.pos[[p]]@positionVector$start))
    }
    tf <- rep(x@bestMatch[[p]]@aligns[[valid]]@TF@name , length(vector.pos[[p]]@positionVector$start))
    chr <- as.character(vector.pos[[p]]@positionVector$chr)
    seq.start <- vector.pos[[p]]@positionVector$start
    seq.end <- vector.pos[[p]]@positionVector$end
    seq.strand <- vector.pos[[p]]@positionVector$motiv.strand
    chr.start <- NULL
    chr.end <- NULL
    chr.strand <- NULL
    #get start, end and chr
    for(j in 1:length(seq.start))
    {
      chr.start <- c(chr.start, y@motifList[[i]]@alignList[[j]]@start)
      chr.end <- c(chr.end, y@motifList[[i]]@alignList[[j]]@end)
      chr.strand <- c(chr.strand, y@motifList[[i]]@alignList[[j]]@strand)
    }
    chr.length = chr.end[1]-chr.start[1]
    start <- NULL
    end <- NULL
    #compute real genomic position
    for (k in 1:length(seq.strand))
    {
      start.tmp <- chr.start[k] + seq.start[k]
      end.tmp <- chr.end[k] - chr.length + seq.end[k]
      start <- c(start, start.tmp)
      end <- c(end, end.tmp)
    }
    ranges <- IRanges(start, end)
    data <- rbind( data, RangedData( ranges, chr, strand=chr.strand, name=mot, motif=tf) )
    p=p+1
  }
  return(data)
}
