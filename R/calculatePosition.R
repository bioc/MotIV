
#####CORRECTION#####
alignmentCorrection <- function (sequence, match, gadem.strand, motiv.strand)
{	
	seq1 <- strsplit(sequence,"")[[1]]
	seq2 <- strsplit(match,"")[[1]]	
	n <- c(0,0)
	for (i in 1:2)
	{
		n1=countGap(seq1, "-")
		n2=countGap(seq2, c("N","-"))
		n[i]=n2-n1
		seq1 <- rev(seq1)
		seq2 <- rev(seq2)
	}
	if (gadem.strand=="+" && motiv.strand=="+") 
	{
		start=n[1]-1
		end=length(seq1[seq1!="-"])-n[2]-2
	}  else if (gadem.strand=="+" && motiv.strand=="-")   { 
		start=n[1]-1
		end=length(seq1[seq1!="-"])-n[2]-2
	}   else if (gadem.strand=="-" && motiv.strand=="+")   { 
		start=-length(seq1[seq1!="-"])+n[2]+1
		end=-n[1]
	}   else if (gadem.strand=="-" && motiv.strand=="-")   { 
		start=-length(seq1[seq1!="-"])+n[2]
		end=-n[1]-1
	} else {
		start=0
		end=0}
	return(c(start, end))
}

######POSITION######
calculatePositionVector <- function (motiv, gadem, group, correction)
{	
	seq.length <- gadem@motifList[[1]]@alignList[[1]]@end-gadem@motifList[[1]]@alignList[[1]]@start
	motiv.count=1;pos.count=1
	pos <- list()
	motivnames <- names(motiv@input)
	
	for (j in 1:length(gadem@motifList))
	{
		sim.pos <- NULL
		start <- NULL
		end <- NULL
		gadem.strand <- NULL
		motiv.strand <- NULL
		seq <- NULL
		chr <- NULL
		
		
		if (gadem@motifList[[j]]@name %in% motivnames)
		{
			current.motiv=which(names(motiv)==gadem@motifList[[j]]@name)
			valid.motiv <- motiv@bestMatch[[current.motiv]]@valid
			if (correction)
			{
				seq.motif <- motiv@bestMatch[[current.motiv]]@aligns[[valid.motiv]]@sequence
				match.motif <- motiv@bestMatch[[current.motiv]]@aligns[[valid.motiv]]@match
				motiv.strand.tmp <- motiv@bestMatch[[current.motiv]]@aligns[[valid.motiv]]@strand
			}
			list.pwm <- list(gadem@motifList[[j]]@pwm)
			names(list.pwm) <- gadem@motifList[[j]]@name
			for (i in 1:length(gadem@motifList[[j]]@alignList))
			{
				position <- gadem@motifList[[j]]@alignList[[i]]@pos
				gadem.strand.tmp <- gadem@motifList[[j]]@alignList[[i]]@strand
				if (correction)
				{
					correctionFactor=alignmentCorrection(seq.motif, match.motif, gadem.strand.tmp, motiv.strand.tmp)
				}  else  {
					correctionFactor=c(0,0)
				}
				start <- c(start, position + correctionFactor[1])	
				end <- c(end,  position + correctionFactor[2])
				
				gadem.strand <- c(gadem.strand, gadem.strand.tmp)
				motiv.strand <- c(motiv.strand, motiv@bestMatch[[current.motiv]]@aligns[[valid.motiv]]@strand )
				seq <- c(seq, gadem@motifList[[j]]@alignList[[i]]@seqID)
				chr <- c(chr, gadem@motifList[[j]]@alignList[[i]]@chr)
			}
			vector.pos <- data.frame(chr=chr,  start=start, end=end, gadem.strand=gadem.strand, motiv.strand=motiv.strand, seq=seq)
			
			if (j>1 && any(similarity(pos) %in% motiv@bestMatch[[current.motiv]]@similarity) && group)
			{
				whichsimilarity=which(similarity(pos)%in% motiv@bestMatch[[current.motiv]]@similarity)
				pos[[whichsimilarity]] <- new("position", motifName=c(pos[[whichsimilarity]]@motifName, gadem@motifList[[j]]@name), positionVector=rbind(pos[[whichsimilarity]]@positionVector, vector.pos), pwm=c(pos[[whichsimilarity]]@pwm, list.pwm), similarity=pos[[whichsimilarity]]@similarity )
			} else   {
				pos[[pos.count]] <- new("position", motifName=gadem@motifList[[j]]@name, positionVector=vector.pos, pwm=list.pwm, similarity=motiv@bestMatch[[current.motiv]]@similarity )
				pos.count=pos.count+1
			}
		}
	}
	return(pos)
}

######################
countGap <- function (seq, char)
{	
	i=1
	while (seq[i] %in% char)
	{  i=i+1 }
	return(i-1)
}