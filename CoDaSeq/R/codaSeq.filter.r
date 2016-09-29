# @title Filter compositional dataset for 0 values and abundance
#
# @ description
# \code{codaSeq.filter} returns a reduced able of counts where the samples must contain a minimum number of reads, and OTUs must be found with a minimum abundance in all remaining samples.
#
# \emph{Note:} filters by min.prop first, min.fraction second, min.prop third.
# Requires numeric data only.
# \emph{Note:} There is a parallel function (to be implemented) that
# reduces the metadata with the same filtering parameters
#
# @param x A matrix or dataframe containing a count table
# @param min.reads The minimum reads per sample. Default=5000
# @param min.prop The minimum proportional abundance of a read in any sample. Default 0.001
# @param max.prop The maximum proportional abundance of a read in any sample. Default 0.025
# @param min.fraction The minimum sample proportion of non-0 reads for each variable
# @param samples.by.row True if rows contain samples, False if rows contain variables
#
# @return A reduced set of data with samples by rows.
#
# @examples
# @ export

codaSeq.filter <- function(x, min.reads=5000, min.prop=0.001, max.prop=1,
  min.occurrence=0, samples.by.row=TRUE){

  if(samples.by.row==FALSE) data <-x
  if(samples.by.row==TRUE) data <- t(x)
  if(length(rownames(data)) == 0) stop("rownames cannot be empty")
  if(length(colnames(data)) == 0) stop("colnames cannot be empty")
  if ( any( round(data) != data ) ) stop("not all values are integers")
  if ( any( data < 0 ) )             stop("one or more values are negative")

  # todo: check for numeric
  data.0 <- data[,which(apply(data,2,sum) > min.reads)]

  d.frac <- apply(data.0, 2, function(x){x/sum(x)})
  data.1 <- data.0[ (which((apply(d.frac,1,max) > min.prop) & (apply(d.frac,1,max) < max.prop))),]
  rm(d.frac)

  data.2 <- data.frame(data.1[which(apply(data.1, 1,
    function(x){length(which(x != 0))/length(x)}) > min.occurrence),])

  return( data.2 )

}
