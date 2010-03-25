exportAsTransfacFile <- function (x, file)
{ 
  if (class(x)!="motiv")
  {stop("motiv must be an object of class motiv.")}

  transfac <- file (paste(file, "_matched.txt", sep=""), "w") #PWMs found
  pairs <-   file (paste(file, "_match_pairs.txt", sep=""), "w") #alignments

  for (i in 1:length(x@input))
  {
    writeChar(paste(">\t", names(x@input)[i], "\n", sep=""), pairs, eos=NULL)
    for (j in 1:length(x@bestMatch[[i]]@aligns))
    {
      writeChar(paste(x@bestMatch[[i]]@aligns[[j]]@TF@name, "\t", format(x@bestMatch[[i]]@aligns[[j]]@evalue, digits=5, scientific=T), "\t", x@bestMatch[[i]]@aligns[[j]]@sequence, "\t", x@bestMatch[[i]]@aligns[[j]]@match, "\n", sep=""), pairs, eos=NULL)
      writeChar(paste("DE\t", x@bestMatch[[i]]@aligns[[j]]@TF@name,"\n", sep=""), transfac, eos=NULL)
      pwm <- format( round(x@bestMatch[[i]]@aligns[[j]]@TF@pwm,4), digits=4)
      for (l in 1:dim(x@bestMatch[[i]]@aligns[[j]]@TF@pwm)[2])
      {
        writeChar(paste(c(l-1, pwm[(4*l-4+1):(4*l)]), c("\t","\t","\t","\t","\n"), sep=""), transfac, eos=NULL)
      }
      writeChar("XX\n", transfac, eos=NULL)
    }

  }
  close(transfac)
  close(pairs)
  print(paste(file,"_matched_transfac.txt created.", sep=""))
  print(paste(file,"_match_pairs.txt created.", sep=""))
}