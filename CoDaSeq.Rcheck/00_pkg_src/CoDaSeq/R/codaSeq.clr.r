#####################
# Center log-ratio function
# returns Center log-ratio transform of data
# @param x A matrix or dataframe with samples by row or column
# @param samples.by.row TRUE if samples by row, FALSE if samples by column
# @return Center log-ratio transform of the data
#####

codaSeq.clr <- function(x, samples.by.row=TRUE){
  if(min(x) == 0) stop("0 values must be replaced, estimated, or eliminated")
  if(samples.by.row == TRUE) margin=1
  if(samples.by.row == FALSE) margin=2

  return ( t(apply(x, margin, function(x){log(x) - mean(log(x))})) )
}
