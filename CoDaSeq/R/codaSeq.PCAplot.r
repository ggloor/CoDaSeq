codaSeq.PCAplot <- function(pcx, plot.groups=FALSE, plot.circles=FALSE, plot.loadings=TRUE,
  grp.col=NULL, grp=NULL, loadings.grp=NULL, loadings.col=NULL, main=""){

    if((attributes(pcx)$class == "prcomp") == FALSE) stop("please use the prcomp function for the SVD")


    if(plot.groups == FALSE) {
      warning("will not plot circles")
      plot.circles=FALSE
      grp.col="black"
      grp=list(c(1:nrow(pcx$x)))
    }

    if(plot.groups == TRUE) {
      plot.circles=plot.circles
      if(is.null(grp.col)) stop("please supply a vector of group colours")
      if(is.null(grp)) stop("please supply a list of sample groupings")
      if(is.list(grp) == FALSE) stop("please supply a list of sample groupings")
      grp.col=grp.col
      grp=grp
    }

         plot.title.l = main
         plot.title.s = main

	mvar <- sum(pcx$sdev^2)
	PC1 <- paste("PC1: ", round(sum(pcx$sdev[1]^2)/mvar, 3))
	PC2 <- paste("PC2: ", round(sum(pcx$sdev[2]^2)/mvar, 3))

	plot(pcx$x[,1], pcx$x[,2], pch=19, xlab=PC1, ylab=PC2,
		col=rgb(0,0,0,0.1),  main= plot.title.s,
		xlim=c(max(abs(pcx$x[,1])) * -1.2, max(abs(pcx$x[,1])) * 1.2),
		ylim=c(max(abs(pcx$x[,2])) * -1.2, max(abs(pcx$x[,2])) * 1.2)
		)
	abline(h=0, lty=2, lwd=2, col=rgb(0,0,0,0.3))
	abline(v=0, lty=2, lwd=2, col=rgb(0,0,0,0.3))

    if(plot.loadings==TRUE){
    if( is.null(loadings.col)){
       scale.factor.1 <- max(abs(pcx$x[,1]))/max(abs(pcx$rotation[,1]))
       scale.factor.2 <- max(abs(pcx$x[,2]))/max(abs(pcx$rotation[,2]))
 	   points(pcx$rotation[,1]*scale.factor.1, pcx$rotation[,2]*scale.factor.2, pch=19,
		col=rgb(0,0,0,0.05), cex=0.5)
    }
    if (!is.null(loadings.col)){
        if(is.null(loadings.col)) stop("please provide a vector of loadings colours")
        if(length(loadings.col) != length(loadings.grp)) stop("number of loadings groups and number of colours must match")
        scale.factor.1 <- max(abs(pcx$x[,1]))/max(abs(pcx$rotation[,1]))
        scale.factor.2 <- max(abs(pcx$x[,2]))/max(abs(pcx$rotation[,2]))
        for(j in 1:length(loadings.grp)){
 	        points(pcx$rotation[,1][loadings.grp[[j]]]*scale.factor.1, pcx$rotation[,2][loadings.grp[[j]]]*scale.factor.2, pch=19, cex=0.5, col=loadings.col[j])

        }
    }
    }



	if(plot.circles == TRUE){
		grp.col=grp.col
		grp = grp
		for(i in 1:length(grp)){
		if (length(pcx$x[grp[[i]],1])  > 2) {
		dataEllipse(pcx$x[grp[[i]],1], pcx$x[grp[[i]],2],
			levels=c(0.67), center.cex=FALSE, plot.points=FALSE, add=TRUE, col=grp.col[i],
			fill = TRUE, fill.alpha = 0.1, pch=19,
			xlim=c(min(pcx$x[,1]) *1.2, max(pcx$x[,1])) *1.2,
			ylim=c(min(pcx$x[,2]) *1.2, max(pcx$x[,2])) *1.2
			)
		}
		graphics::text(pcx$x[grp[[i]],1], pcx$x[grp[[i]],2], labels=rownames(pcx$x)[grp[[i]]],
           col=grp.col[i])

		}
	} else if(plot.circles == FALSE){
		grp.col=grp.col
		grp = grp
		for(i in 1:length(grp)){
		graphics::text(pcx$x[grp[[i]],1], pcx$x[grp[[i]],2], labels=rownames(pcx$x)[grp[[i]]],
		col=grp.col[i])
	    }
	}
}

