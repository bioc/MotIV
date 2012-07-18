
motifMatch <- function (inputPWM, database=jaspar, DBscores=jaspar.scores, cc="PCC", align="SWU", top=5, go=1, ge=0.5) {
	res<-.Call("RmotifMatch", cc=cc, align=align, top=top, go=go, ge=ge, inputPWM=inputPWM, inputDB=database, inputScores=DBscores)
	topx=length(res[[7]])/length(res[[1]])
	if (!is.null(res))
	{
		list1 <- NULL
		for (i in 1:(topx*length(res[[1]])))
		{
			tf <- new("transcriptionFactor", name=res[[2]][i], pwm=res[[3]][[i]])
			sequence.strand <- res[[7]][i]
			match.strand <- res[[8]][i]			
			if (match.strand == "-") 
			{
				sequence <- as.character(reverseComplement(DNAString(res[[5]][i])))
				match <- as.character(reverseComplement(DNAString(res[[6]][i])))
			} else {
				sequence <- res[[5]][i]
				match <- res[[6]][i]
			}			
			alig <- new("alignments", TF=tf, evalue=res[[4]][i], sequence=sequence, match=match, strand=res[[8]][i])	
			list1 <- c(list1, alig)
		}
		
		bestMatch <- NULL
		for (i in 1:length(res[[1]]))
		{
			bestMatch <- c(bestMatch, new("matches", name=res[[1]][[i]], aligns=list1[(topx*i-topx+1):(topx*i)], similarity=res[[1]][[i]], valid=1))
		}
		
		motiv <- new("motiv", input=inputPWM, bestMatch=bestMatch, argv=c(cc=cc, align=align, top=top, go=go, ge=ge))
		return(motiv)
	}
}


