generateDBScores <- function (inputDB, cc="PCC", align="SWU", nRand=1000, go=1, ge=0.5)
{
  res <- .Call("RgenerateScoresDB", cc=cc, align=align, go=go, ge=ge, nRand=nRand, inputPWM=inputDB)
  colnames(res) <- c("length1", "length2", "mean", "var", "count", "min", "max")
  return(res)
}

writeDBScores <- function (x, file)
{
  write.table(x=x, file=file, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
}

readDBScores <- function (file)
{
  scores <- as.matrix(read.csv(file=file, sep="\t", header=FALSE))
  return(scores)
}
