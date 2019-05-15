
#####CORRECTION#####

alignmentCorrection <- function (sequence, match)
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
	correction.table=matrix(ncol=5, nrow=2)
	correction.table[,1]=c(start=n[1]-1, end=length(seq1[seq1!="-"])-n[2]-2)
	correction.table[,2]=c(start=n[1]-1, end=length(seq1[seq1!="-"])-n[2]-2)
	correction.table[,3]=c(start=-length(seq1[seq1!="-"])+n[2]+1, end=-n[1])
	correction.table[,4]=c(start=-length(seq1[seq1!="-"])+n[2], end=-n[1]-1) #correction of -1bp to match with getSeq
	correction.table[,5]=c(0,0)
	return(correction.table)
}



######POSITION######
calculatePositionVector <- function (motiv, gadem, group, correction)
{	
	motiv.count=1;pos.count=1
	pos <- list()
	sequencesLength=sapply(gadem@motifList, function(x){sapply(x@alignList, function(x){x@end})-sapply(x@alignList, function(x){x@start})})
		
	for (j in seq(1,length(motiv)))
	{
		
			current.motiv=which( names(gadem)==names(motiv)[j])

			#current.motiv=which(names(motiv)==gadem@motifList[[j]]@name)
			valid.motiv <- motiv@bestMatch[[j]]@valid
			PWMs<- list(gadem@motifList[[current.motiv]]@pwm)
			names(PWMs) <- gadem@motifList[[current.motiv]]@name
							
			seq.motif <- motiv@bestMatch[[j]]@aligns[[valid.motiv]]@sequence
			match.motif <- motiv@bestMatch[[j]]@aligns[[valid.motiv]]@match
			gadem.strand <-  sapply(gadem@motifList[[current.motiv]]@alignList, function(x){x@strand}) 
			motiv.strand <- motiv@bestMatch[[j]]@aligns[[valid.motiv]]@strand

			#alignments <- data.frame(sequence=seq.motif, match=match.motif, strand1=gadem.strand, strand2=motiv.strand)
			
			chr <- sapply(gadem@motifList[[current.motiv]]@alignList, function(x){x@chr})  #gsub("chr","", sapply(gadem@motifList[[current.motiv]]@alignList, function(x){x@chr}) )
			position <- sapply(gadem@motifList[[current.motiv]]@alignList, function(x){x@pos + x@start })
			
			 # peakPos <- sapply(gadem@motifList[[current.motiv]]@alignList, function(x){x@pos })
			
			peakPos <- (-sequencesLength[[current.motiv]])/2 + sapply(gadem@motifList[[current.motiv]]@alignList, function(x){x@pos })
			
			seqID <- sapply(gadem@motifList[[current.motiv]]@alignList, function(x){x@seqID})
			
			lengthPeak <- sapply(gadem@motifList[[current.motiv]]@alignList, function(x){x@end-x@start})
			
			if(correction)
			{
				correction.table <- alignmentCorrection(as.character(seq.motif), as.character(match.motif))
				correctionFactor=correction.table[,unlist(sapply(paste(gadem.strand,motiv.strand, sep=""), function(x){switch(x, "++"=1, "+-"=2, "-+"=3, "--"=4, 5)}))] 
				
				data.ir <- IRanges(start=position + correctionFactor[1,], end=position + correctionFactor[2,])
				data.gr <- GRanges(chr, data.ir, gadem.strand, seqID=seqID, lengthPeak= lengthPeak, peakPos=  (2*peakPos+ correctionFactor[1,] + correctionFactor[2,] )/2) 
			} else {
				data.ir <- IRanges(start=position + sapply(gadem@motifList[[current.motiv]]@alignList, function(x){x@start}), end=position + sapply(gadem@motifList[[current.motiv]]@alignList, function(x){x@start}) + length(gadem@motifList[[current.motiv]]@consensus ))
				data.gr <- GRanges(chr, data.ir, gadem.strand,  seqID=seqID, lengthPeak=lengthPeak, peakPos=peakPos)
			}

			if (pos.count>1 && any(similarity(pos) %in% motiv@bestMatch[[j]]@similarity) && group)
			{
				whichsimilarity=which(similarity(pos)%in% motiv@bestMatch[[j]]@similarity)
				pos[[whichsimilarity]] <- new("position", motifName=c(pos[[whichsimilarity]]@motifName, gadem@motifList[[current.motiv]]@name), positionVector=c(pos[[whichsimilarity]]@positionVector, data.gr), pwm=c(pos[[whichsimilarity]]@pwm, PWMs), similarity=pos[[whichsimilarity]]@similarity )
			} else   {
				pos[[pos.count]] <- new("position", motifName=gadem@motifList[[current.motiv]]@name, positionVector=data.gr, pwm=PWMs, similarity=motiv@bestMatch[[j]]@similarity)
				pos.count=pos.count+1
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
