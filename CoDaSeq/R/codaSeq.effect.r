#' Estimate effect size between two distributions
#'
#' Calculates relatively robust measures of standardized mean difference (effect) and dispersion from a matrix or dataframe of numbers. The output value is smaller thant Cohen's d by a factor of 1.418 when comparing two Normal distributions. There is an option to normalize the effect size by that factor.
#' @param x A numerical matrix with samples by column
#' @param conds A vector denoting group membership
#' @param corrected Whether to scale to Cohen's d or not, default is FALSE
#' @return returns a vector of effect sizes
#' @export
#' @examples
#' # make a synthetic dataset
#' d <- c(rnorm(100,0,1), rnorm(100,2,1))
#' e <- c(rnorm(100,2,1), rnorm(100,0,1))
#'  de <- rbind(d,e)
#' conds <- c(rep("A", 100), rep("B",100))
#' # values should be approximately -2 and 2
#' codaSeq.effect(de, conds, corrected=TRUE)


codaSeq.effect <- function(x, conds,corrected=FALSE){
    conds <- as.factor(conds)
    levels <- levels(conds)
    levels <- vector("list", length(levels))
    names(levels) <- levels(conds)
    for ( l in levels( conds ) ) {
        levels[[l]] <- which( conds == l )
        if ( length( levels[[l]] ) < 2 ) stop("condition level '",l,"' has less than two replicates")
    }
    return.data <- vector()
    for(i in 1:nrow(x)){
        return.data[i] <- (private.dnef(as.numeric(x[i,levels[[1]] ]), as.numeric(x[i,levels[[2]] ])))
    }
    if(corrected == FALSE) return(return.data)
    if(corrected == TRUE) return(return.data*1.418)
}

# vectorized maximum absolute deviation
private.mxad <- function(x,y){
  pmax(abs(x - sample(x, size=10*length(x), replace=T)), abs(y - sample(y,size=10*length(y), replace=T)))
}

private.dnef <- function(a,b){
   return(median( (a - b)/ (private.mxad(a,b)) ) )
}

