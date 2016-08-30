#####################
# Center log-ratio function
# returns Center log-ratio transform of data
# @param x A matrix or dataframe with samples by row or column
# @param samples.by.row TRUE if samples by row, FALSE if samples by column
# @return Center log-ratio transform of the data with samples by row
#####

codaSeq.clr <- function(x, IQLR=FALSE, samples.by.row=TRUE){
  if(min(x) <= 0) stop("only positive real values permitted")
  if(samples.by.row == TRUE) margin=1
  if(samples.by.row == FALSE) margin=2

  if (IQLR == FALSE){
      return( t(apply(x, margin, function(x){log(x) - mean(log(x))})) )
  } else if (IQLR == TRUE){
      reads.clr <- t(apply(x, margin, function(x){log(x) - mean(log(x))}))
      reads.var <- apply(reads.clr, 2, var)
      reads.qtl <- quantile(unlist(reads.var))
      mid.set <- which(reads.var < (reads.qtl[4]) & reads.var > (reads.qtl[2]))
      return(t(apply(x, margin, function(x) log(x) - mean(log(x[mid.set])))))
  }
}
