draw.grey.boxes <- function(at) {
   xlim <- par('usr')[1:2] # gets the limit of the bounding box
   for ( a in 1:length(at) ) {
       if ( a %% 2 == 1 ) next
       rect( xlim[1], a - 0.5, xlim[2], a + 0.5, border=NA, col=rgb(0,0,0,0.1) )
   }

}

codaSeq.stripchart <- function(
    aldex.out=NULL, group.table=NULL, group.label=NULL, sig.method="we.eBH",
    x.axis="effect", sig.cutoff=0.05, cex=0.8, main=NULL, mar=c(4,12,4,0.5),
    do.ylab=TRUE, heir=FALSE, heir.base=NULL)
  {
  # aldex.out is the data frame of observations to be plotted
  # group.table taxon information. Each column is a taxon name or function at a specific level
  # group.label is the column name from the taxon or function file
  # x.axis should be either diff.btw (difference) or effect
  # cex controls plot size
  # p.cutoff is the BH corrected wilcoxon rank test statistic
  # effect.cutoff is the BH corrected wilcoxon rank test statistic
  # do.ylab controls whether to add y axis label to plot for pseudo-time series display
  # copy the aldex output data frame

  if(is.null(aldex.out)) stop("please supply an appropriate input from ALDEx2")
  if(is.null(group.table)) stop("please supply an appropriate group table")
  if(is.null(group.label)) stop("please supply an appropriate group label")
  if(sig.method != "effect" & sig.cutoff > 0.1) stop("p value cutoff not realistic")

  non.sig <- list() # not sig
  sig <- list() # significant by effect or by p

  if(heir == FALSE){
    aldex.out <- data.frame(aldex.out, group.table[rownames(aldex.out), group.label])
    colnames(aldex.out)[ncol(aldex.out)] <- group.label
    aldex.out <- droplevels(aldex.out)

    # make a list of unique groups in the dataset
    # there may be empty groups.set, ignore these for now
    groups.set <- sort(unique(aldex.out[[group.label]]), decreasing=TRUE)

	for(i in 1:length(groups.set)){
		grp.nms <- rownames(aldex.out)[aldex.out[,group.label] == groups.set[i]]

		if(sig.method == "effect"){
		  non.sig[[as.character(groups.set[i])]] <- aldex.out[grp.nms,x.axis][abs(aldex.out[grp.nms, sig.method]) <= sig.cutoff]
		  sig[[as.character(groups.set[i])]] <- aldex.out[grp.nms,x.axis][abs(aldex.out[grp.nms, sig.method]) > sig.cutoff]
		}else{
		  non.sig[[as.character(groups.set[i])]] <- aldex.out[grp.nms,x.axis][abs(aldex.out[grp.nms, sig.method]) >= sig.cutoff]
		  sig[[as.character(groups.set[i])]] <- aldex.out[grp.nms,x.axis][abs(aldex.out[grp.nms, sig.method]) < sig.cutoff]
		}
	}


  }else if(heir == TRUE){
  	if(is.null(heir.base)) stop("please supply the index or name of the heirarchy base")
    nms <- group.table[group.table[,heir.base] %in% rownames(aldex.out),]
    groups.set <- sort(unique(nms[[group.label]]), decreasing=TRUE)

	for(i in 1:length(groups.set)){
    	grp.nms <- nms[,heir.base][nms[,group.label] == groups.set[i]]

	if(sig.method == "effect"){
		  non.sig[[as.character(groups.set[i])]] <- aldex.out[grp.nms,x.axis][abs(aldex.out[grp.nms, sig.method]) <= sig.cutoff]
		  sig[[as.character(groups.set[i])]] <- aldex.out[grp.nms,x.axis][abs(aldex.out[grp.nms, sig.method]) > sig.cutoff]
		}else{
		  non.sig[[as.character(groups.set[i])]] <- aldex.out[grp.nms,x.axis][abs(aldex.out[grp.nms, sig.method]) >= sig.cutoff]
		  sig[[as.character(groups.set[i])]] <- aldex.out[grp.nms,x.axis][abs(aldex.out[grp.nms, sig.method]) < sig.cutoff]
		}
    }
  }

  # generate a y axis plotting limit a bit larger than needed
  ylim<-c(length(groups.set) - (length(groups.set)+0.5), length(groups.set) + 0.5)

  x.axis <- x.axis # we find the effect or diff.btw columns to be most useful
  xlim = c(min(-1 * max(abs(aldex.out[,x.axis]))), max(aldex.out[,x.axis]))

  # basically can call the different significant groups.set within strip chart
  par(mar=mar, las=1, cex=cex)
if(do.ylab == TRUE) { stripchart(non.sig,
    col=c(rgb(0,0,0,0.2),rgb(0,0,0,0.3)), method="jitter", pch=19, xlim=xlim, xlab=x.axis, main=main)
    }
if(do.ylab == FALSE) {stripchart(non.sig,
    col=rgb(0,0,0,0.3), method="jitter", pch=19, xlim=xlim, xlab=x.axis, main=main, yaxt="n")
    }

draw.grey.boxes(as.vector(groups.set))
sig.cex=cex+0.2

stripchart(sig, col=rgb(1,0,0,0.3),method="jitter", pch=19, add=T, cex=sig.cex)

abline(v=0, lty=2, col=rgb(0,0,0,0.2),lwd=2)
abline(v=1, lty=3, col=rgb(0,0,0,0.2),lwd=2)
abline(v= -1, lty=3, col=rgb(0,0,0,0.2),lwd=2)

}
