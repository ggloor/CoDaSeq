#' Identifying sample outliers
#'
#' Returns a list of proportional contibution to group variance,
#' sample names that are outliers, and sample names that are not outliers.
#'
#' @param x A matrix or data frame with CLR-transformed values, 
#'   with samples by row; ideally from `codaSeq.clr()` or `aldex.clr()`.
#' @param plot.me A logical value determining if a histogram of the variance 
#'   contribution per sample should be plotted (default: `TRUE`).
#' @param col RGB values for your selection of colours (default: `rgb(1,0,0,0.3)`).
#'
#' @return Returns a list with three elements: `sample.var`, containing 
#'   proportional variance contributions for each sample, `bad`, containing
#'   outlier samples, and `good`, containing non-outlier samples.
#' 
#' @export
#'
#' @details Samples must be grouped, but this approach makes no sense when used across 
#'   groups. If you do not know if you have natural groups, ignore this step and
#' 	 examine your data by PCA. Outliers are defined as those contributing 
#' 	 more than the median plus twice the interquartile range of the sample 
#' 	 variance to the total.
#'
#' @examples
#' # load example HMP data from CoDaSeq package
#' data("ak_op")  # feature table: 4347 OTU x 30 samples (15x ak, 15x op)
#' 
#' # Bayesian-multiplicative replacement of count-zeros 
#' clr.input<-cmultRepl(t(ak_op),label = "0",
#'                      method = "CZM",output = "p-counts")
#' 
#' # CLR-transformation using codaSeq.clr
#' clr.data<-codaSeq.clr(t(clr.input), IQLR = FALSE,
#'                       aitch = FALSE, samples.by.row = TRUE)
#' # get outliers for keratinised gingiva samples
#' outliers.ak<-codaSeq.outlier(clr.data[1:15,])
#'
#' # get outliers for oral plaque samples
#' outliers.ak<-codaSeq.outlier(clr.data[16:30,])
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
