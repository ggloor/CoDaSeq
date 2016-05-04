# @title Filter compositional dataset for 0 values and abundance
#
# @ description
# \code{codaSeq.filter} returns a reduced able of counts where the
# samples must contain a minimum number of reads, and OTUs must be
# found with a minimum abundance in all remaining samples.
# \emph{Note:} filters by min.prop first, min.fraction second, min.prop third.
# Requires numeric data only.
# \emph{Note:} There is a parallel function (to be implemented) that
# reduces the metadata with the same filtering parameters
#
# @param x A matrix or dataframe containing a count table
# @param min.reads The minimum reads per sample. Default=5000
# @param min.prop The minimum proportional abundance of a read in any sample. Default 0.001
# @param min.fraction The minimum sample proportion of non-0 reads for each variable
# @param samples.by.row True if rows contain samples, False if rows contain variables
#
# @return A reduced set of data with samples by rows
#
# @examples
# to do add taxonomy filter
# @ export
codaSeq.filter <- function(x, y=tax.vector, min.reads=5000, min.prop=0.001,
  min.occurrence=0, samples.by.row=TRUE){
  if(samples.by.row==FALSE) data <-x
  if(samples.by.row==TRUE) data <- t(x)
  # todo: check for numeric
  data.0 <- data[,which(apply(data,2,sum) > min.reads)]

  d.frac <- apply(data.0, 2, function(x){x/sum(x)})
  data.1 <- data.0[which(apply(d.frac, 1, max) > min.prop),]
  rm(d.frac)

  data.2 <- data.frame(data.1[which(apply(data.1, 1,
    function(x){length(which(x != 0))/length(x)}) > min.occurrence),])

  return( data.2 )
}

# @title Center log-ratio function
#
# @description
# \code{codaSeq.clr} returns a matrix of center log-ratio transformed data
# with samples by row
# equivalent to log(x/gx) for every value where gx is the geomtric mean of the vector X
# \emph{Note:} natural log is used for biplots and other exploratory analyses
# @param x A matrix or dataframe with samples by row or column
# @param samples.by.row TRUE if samples by row, FALSE if samples by column
# @return Center log-ratio transform of the data
codaSeq.clr <- function(x, samples.by.row=TRUE){
  if(min(x) == 0) stop("0 values must be replaced, estimated, or eliminated")
  if(samples.by.row == TRUE) margin=1
  if(samples.by.row == FALSE) margin=2

  return ( t(apply(x, margin, function(x){log(x) - mean(log(x))})) )
}

# @title Identifying sample outliers
#
# @description
# \code{codaSeq.outlier} returns a list of proportional contribution to group variance,
# sample names that are outliers, and sample names that are not outliers
# \emph{Note:} Samples must be grouped. This approach makes no sense across
# groups. If you do not know if you have natural groups, ignore this step and
# examine your data by PCA
# outliers are defined as those contributing greater than the median plus twice the
# interquartile range of the sample variance to the total
# @param x A matrix or dataframe with clr transformed values, with samples by row
# @param plot.me A logical value determining if a histogram should be plotted of the
# variance contribution per sample
# @return list
# @sample.var: proportional variance contributions for each sample
# @bad: rownames of outlier samples
# @good: rownames of non-outlier samples
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

# @title Rarefy samples to common sample depth

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



