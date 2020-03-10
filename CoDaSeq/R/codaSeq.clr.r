#' Center log-ratio function
#'
#' Returns a matrix of center log-ratio transformed data with samples
#'   by row. Each value is equivalent to log(x/gx) where gx is
#'   the geometic mean of the row vector X.
#' @param x A matrix or dataframe with samples by row or column
#' @param IQLR the geometric mean computed on all features if FALSE,
#'   or on the set of features with variance between the first and
#'   third quartile if TRUE. To be used when the data is not centered.
#' @param aitch returns the values drawn from a digamma distribution using
#'   an uniformative prior. Approximates the expected value of the clr
#'   values and is much faster than using aldex.effect
#' @param samples.by.row TRUE if samples by row, FALSE if samples by column
#' @return returns a matrix of clr tranformed values with samples
#'   in the rows and variables in columns
#' @details the natural logarithm is used.
#' @author Greg Gloor, Andrew Fernandes, Jean Macklaim, Wallace Chan
#' @seealso \code{\link{codaSeq.filter}},
#' \code{\link{codaSeq.outlier}},
#' \code{\link{codaSeq.rarefy}},
#' \code{\link{codaSeq.phi}}
#'
#' @examples
#' # get dataset from ALDEx2 package
#' data(selex)
#' # convert to clr with an uninformative prior
#' # use only the variables with mid-quartile variance as denominator
#' # output will have samples by rows
#' selex.clr <- codaSeq.clr(selex + 0.5, IQLR=TRUE, aitch=TRUE, samples.by.row=FALSE)

#####

codaSeq.clr <- function(x, IQLR=FALSE, aitch=FALSE, samples.by.row=TRUE){
  if(min(x) < 0) stop("only positive real values permitted")
  if (!is.vector(x[,1], mode="numeric") ) stop("counts must be supplied as numbers")
  if ( any( x < 0 ) ) stop("counts cannot be negative")
  if(samples.by.row == TRUE) margin=1
  if(samples.by.row == FALSE) margin=2

  if (aitch == FALSE){
  if (IQLR == FALSE){
      if(samples.by.row == T) return( t(apply(x, margin, function(x){log(x) - mean(log(x))})) )
      if(samples.by.row == F) return( apply(x, margin, function(x){log(x) - mean(log(x))}) )
  } else if (IQLR == TRUE){
      reads.clr <- t(apply(x, margin, function(x){log(x) - mean(log(x))}))
      reads.var <- apply(reads.clr, 2, var)
      reads.qtl <- quantile(unlist(reads.var))
      mid.set <- which(reads.var < (reads.qtl[4]) & reads.var > (reads.qtl[2]))
      if(samples.by.row == F) return(apply(x, margin, function(x) log(x) - mean(log(x[mid.set]))))

      if(samples.by.row == T) return(t(apply(x, margin, function(x) log(x) - mean(log(x[mid.set])))))
  }
  }
  if (aitch == TRUE){
    aitchison.mean <- function( n, log=TRUE ) {

        # Input is a vector of non-negative integer counts.
        # Output is a probability vector of expected frequencies.
        # If log-frequencies are requested, the uninformative subspace is removed.

        a <- n + 0.5
        sa <- sum(a)

        log.p <- digamma(a) - digamma(sa)
        log.p <- log.p - mean(log.p)

        if ( log ) return(log.p)

        p <- exp( log.p - max(log.p) )
        p <- p / sum(p)
        return(p)
    }

    if(samples.by.row == FALSE){
        return(apply(x, margin, aitchison.mean))
    }
    if(samples.by.row == TRUE) x <- t(x)
        return(apply(x, margin, aitchison.mean))
    }
}
