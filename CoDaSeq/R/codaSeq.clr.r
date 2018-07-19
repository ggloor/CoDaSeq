#' Center log-ratio function
#'
#' Returns a matrix of center log-ratio transformed data with samples
#'   by row. Each value is equivalent to log(x/gx) where gx is
#'   the geometic mean of the row vector X.
#' @param x A matrix or dataframe with samples by row or column
#' @param IQLR the geometric mean computed on all features if FALSE,
#'   or on the set of features with variance between the first and
#'   third quartile if TRUE. To be used when the data is not centered.
#' @param samples.by.row TRUE if samples by row, FALSE if samples by column
#' @return returns a matrix of clr tranformed values with samples
#'   in the rows and variables in columns
#' @details the natural logarithm is used.
#' @author Greg Gloor, Jean Macklaim, Wallace Chan
#' @seealso \code{\link{codaSeq.filter}},
#' \code{\link{codaSeq.outlier}},
#' \code{\link{codaSeq.rarefy}},
#' \code{\link{codaSeq.phi}}
#'
#' @examples
#' # get dataset from ALDEx2 package
#' data(selex)
#' # convert to clr with an uninformative prior
#' use only the variables with mid-quartile variance as denominator
#' output will have samples by rows
#' selex.clr <- codaSeq.clr(selex + 0.5, IQLR=TRUE, samples.by.row=FALSE)

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
