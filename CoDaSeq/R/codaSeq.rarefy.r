#' Subsample a frequency table without replacement
#'
#' Returns a table of counts where samples have been sampled without replacement.
#' This is included for compatibility, but in general is a bad idea, since it
#' results in a loss of information, and distorts the underlying data somewhat.
#'
#' @param x A matrix or dataframe containing a count table.
#' @param n The desired target number of reads per sample (Default: `1000`).
#' @param samples.by.row A logical value indicating if samples are in rows (`TRUE`)
#'   or in columns (`FALSE`) (Default: `TRUE`).
#'
#' @return Returns a matrix in the same orientation as the original with counts
#'  per OTU reduced to the common sampling depth. This is a constant sum operation.
#'  
#' @export
#' 
#' @author Greg Gloor, Jean Macklaim, Wallace Chan
#'
#' @examples
#' #' # load example HMP data from CoDaSeq package
#' data("ak_op")  # feature table: 4347 OTU x 30 samples (15x ak, 15x op)
#' 
#' # calculate total number of reads in each sample (all are >1000)
#' colSums(ak_op)
#' 
#' # rarefy read table to a depth of 1000 reads
#' codaSeq.rarefy(ak_op, n= 1000, samples.by.row= FALSE)
#' 
#' 
codaSeq.rarefy <- function(x, n=1000, samples.by.row=TRUE){

	rarefy <- function(x,n){
		lx <- length(x)
		vec <- vector()
		rare.vec <- vector()

		for(i in 1:lx){ vec <- c(vec,rep(i, x[i]) ) }
		rare <- sample(vec, n, replace=FALSE)
		for(i in 1:lx){ rare.vec[i] <- length(which(rare == i))}

		return(rare.vec)
	}

    if(samples.by.row == TRUE) margin=1
    if(samples.by.row == FALSE) margin=2
	d.rare <- data.frame(t(apply(x, margin, function(x){rarefy(x,n)})))

    if(samples.by.row == FALSE) {
      rownames(d.rare) <- colnames(x)
      colnames(d.rare) <- rownames(x)
      }
    if(samples.by.row == TRUE) {
      rownames(d.rare) <- rownames(x)
      colnames(d.rare) <- colnames(x)
      }
    if(samples.by.row==TRUE){
      return(d.rare)
      } else return(t(d.rare))
}
