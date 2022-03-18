#####################
# Identifying sample outliers
# returns list of sample outliers
# this is now deprecated because it can be replaced by a 
# much simpler version in progress
#####

codaSeq.outlier <- function(x, plot.me=TRUE, col=rgb(1,0,0,0.3)){

  pcx <- prcomp(x)
  mv <- sum(pcx$sdev^2)

  sample.var <-  apply(pcx$x,1,function(y){sum(y^2/mv)})

  cut <- median(apply(pcx$x,1,function(x){sum(x^2/mv)})) + 2 * IQR(apply(pcx$x,1,function(x){sum(x^2/mv)}))

  bad <- names(which(apply(pcx$x,1,function(x){sum(x^2/mv)}) > cut))
  good <- names(which(apply(pcx$x,1,function(x){sum(x^2/mv)}) <= cut))
  if(plot.me == TRUE){
    hist(sample.var, breaks=100)
    boxplot(sample.var, horizontal=TRUE, col=col, add=TRUE)
    abline(v=cut, lty=2)
  }
  return(list(sample.var=sample.var, bad=bad, good=good) )
}
