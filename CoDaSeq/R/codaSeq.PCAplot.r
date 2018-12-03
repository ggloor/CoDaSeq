codaSeq.PCAplot <- function(pcx, plot.groups=FALSE, plot.circles=FALSE, plot.loadings=TRUE,
  grp.col=NULL, grp=NULL, main=""){
    if((attributes(pcx)$class == "prcomp") == FALSE) stop("please use the prcomp function for the SVD")


    if(plot.groups == FALSE) {
      warning("will not plot circles")
      plot.circles=FALSE
      grp.col="grey"
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

    if(main == NULL){
         plot.title.l = NULL
         plot.title.s = NULL
    } else {
        plot.title.l = paste(main, "\nloadings", sep="")
        plot.title.s = paste(main, "\nsamples", sep="")
    }

	mvar <- sum(pcx$sdev^2)
	PC1 <- paste("PC1: ", round(sum(pcx$sdev[1]^2)/mvar, 3))
	PC2 <- paste("PC2: ", round(sum(pcx$sdev[2]^2)/mvar, 3))

	if(plot.loadings == TRUE){
	  plot(pcx$rotation[,1], pcx$rotation[,2], pch=19, xlab=PC1, ylab=PC2,
		col=rgb(0,0,0,0.1),  main= plot.title.l,
		xlim=c(min(pcx$rotation[,1]) , max(pcx$rotation[,1])) * 1.2,
		ylim=c(min(pcx$rotation[,2]) , max(pcx$rotation[,2])) * 1.2
	  )
	  abline(h=0, lty=2, lwd=2, col=rgb(0,0,0,0.3))
	  abline(v=0, lty=2, lwd=2, col=rgb(0,0,0,0.3))

	}

	plot(pcx$x[,1], pcx$x[,2], pch=19, xlab=PC1, ylab=PC2,
		col=rgb(0,0,0,0.1),  main= plot.title.s,
		xlim=c(min(pcx$x[,1]) *1.2, max(pcx$x[,1])) *1.2,
		ylim=c(min(pcx$x[,2]) *1.2, max(pcx$x[,2])) *1.2
		)
	abline(h=0, lty=2, lwd=2, col=rgb(0,0,0,0.3))
	abline(v=0, lty=2, lwd=2, col=rgb(0,0,0,0.3))

	if(plot.circles == TRUE){
		grp.col=grp.col
		grp = grp
		for(i in 1:length(grp)){
		dataEllipse(pcx$x[grp[[i]],1], pcx$x[grp[[i]],2],
			levels=c(0.67), center.cex=FALSE, plot.points=TRUE, add=TRUE, col=grp.col[i],
			fill = TRUE, fill.alpha = 0.2, pch=19,
			xlim=c(min(pcx$x[,1]) *1.2, max(pcx$x[,1])) *1.2,
			ylim=c(min(pcx$x[,2]) *1.2, max(pcx$x[,2])) *1.2
			)
		}
	} else if(plot.circles == FALSE){
		grp.col=grp.col
		grp = grp
		for(i in 1:length(grp)){
		points(pcx$x[grp[[i]],1], pcx$x[grp[[i]],2], pch=19, col=grp.col[i],
		xlim=c(min(pcx$x[,1]) *1.2, max(pcx$x[,1])) *1.2,
		ylim=c(min(pcx$x[,2]) *1.2, max(pcx$x[,2])) *1.2)
	    }
	}
}

