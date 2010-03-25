
#####TRIM#####

trimPWMedge <- function(x, threshold=1)
{
  res <- list()
  for (i in seq(length(x)))
  {
    for (k in 1:2)
    {
      j=1
      ic  <- x2ic(x[[i]]/sum(x[[i]][,1]))
      if(all(ic < threshold))
      {
        stop("Threshold too high.")
      }
      while (ic[j] < threshold)
      {
        res[[i]] <- x[[i]][,-1]
        j=j+1
      }	
      res[[i]] <- x[[i]][, dim(x[[i]])[2]:1]
    }
  }
  names(res) <- names(x)
  return (res)
}

x2ic <- function(x)
{
  npos <- ncol(x)
  ic <- numeric(length=npos)
  for (i in 1:npos)
  {
    ic[i] <- 2 + sum(sapply(x[, i], function(x) { 
      if (x > 0) { x*log2(x) } else { 0 }
      }))
    }    
    ic
  }
