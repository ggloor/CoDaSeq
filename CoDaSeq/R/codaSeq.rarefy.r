#####################
# Rarefy samples to common sample depth
####

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

    if(samples.by.row == TRUE) margin=2
    if(samples.by.row == FALSE) margin=1
	d.rare <- data.frame(t(apply(x, margin, function(x){rarefy(x,n)})))

    if(samples.by.row == FALSE) {
      rownames(d.rare) <- rownames(x)
      colnames(d.rare) <- colnames(x)
      }
    if(samples.by.row == TRUE) {
      rownames(d.rare) <- colnames(x)
      colnames(d.rare) <- rownames(x)
      }
    return(d.rare)
}
