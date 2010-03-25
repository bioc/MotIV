
readGademPWMFile <- function (file)
{
  table <- read.csv(file=file , fill=T, header=F)
  skiplines <- seq(1, dim(table)[1],5) #skip ">" lines
  table.pwm <- as.character(table[-skiplines, ])
  PWM <- strsplit(table.pwm ,"[\t, ]")

  PWM.names <- as.character(table[skiplines,])
  PWM.names <- sapply(PWM.names,    function(a){strsplit(a,">")[[1]][2]}    , USE.NAMES=F, simplify=T)

  listPWM <- list()
  for (i in seq(length(PWM.names)))
  {
    p=4*(i-1)+1
    listPWM[[i]] <- matrix(as.numeric(c(PWM[[p]][-1], PWM[[p+1]][-1], PWM[[p+2]][-1], PWM[[p+3]][-1])) , nrow=4, byrow=T, dimnames=list(c("A","C","G","T")))
    colnames(listPWM[[i]]) <- 1:(length(listPWM[[i]])/4)
    names(listPWM)[i] <- PWM.names[i]
  }
  return(listPWM)
}

