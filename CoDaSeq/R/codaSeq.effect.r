#####################
# return effect size and MMAD
# requries a vector denoting group membership
#####

# vectorized maximum absolute deviation
codaSeq.mxad <- function(x,y){
  pmax(abs(x - sample(x, size=10*length(x), replace=T)), abs(y - sample(y,size=10*length(y), replace=T)))
}

codaSeq.dnef <- function(a,b){
   return(median( (a - b)/ (codaSeq.mxad(a,b)) ) )
}

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
        return.data[i] <- (codaSeq.dnef(as.numeric(x[i,levels[[1]] ]), as.numeric(x[i,levels[[2]] ])))
    }
    if(corrected == FALSE) return(return.data)
    if(corrected == TRUE) return(return.data*1.418)
}
