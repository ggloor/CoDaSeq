\name{codaSeq.filter}
\alias{codaSeq.filter}
\title{Filter compositional dataset for 0 values and abundance.}
\description{
	Returns a reduced able of counts where the samples must contain a minimum number of
	reads, and OTUs must be found with a minimum abundance in all remaining samples.
	Filters are applied sequentially
}
\usage{
codaSeq.filter <- function(x, min.reads=5000, min.prop=0.001, max.prop=1,
    min.occurrence=0, samples.by.row=TRUE)
}
\arguments{
	\item{x}{
		A matrix or dataframe containing a count table.
	}
	\item{min.reads}{
		The minimum reads per sample. Default=5000.
	}
	\item{min.prop}{
		The minimum proportional abundance of a read in any sample. Default=0.001.
	}
	\item{max.prop}{
		The maximum proportional abundance of a read in any sample. Default=1.
	}
	\item{min.occurrence}{
		The minimum fraction of non-0 reads for each variable in all samples.
	}
	\item{sample.by.row}{
		True if rows contain samples, false if rows contain variables.
	}
}
\details{
	Filters by min.reads first,  min/max.prop second, min.occurrence third.
	Requires numeric data only.
}
\value{
	Returns a dataframe with the following information:
	\item{data.0}{
	}
	\item{data.1}{
	}
	\item{data.2}{
		Returns a reduced vector with filtered samples by rows.
	}
}
\references{
	Please use the citation given by \code{citation(package="CoDaSeq")}
}
\author{
	Greg Gloor, Jean Macklaim, Wallace Chan
}
\seealso{
	\code{\link{codaSeq.clr}},
	\code{\link{codaSeq.outlier}},
	\code{\link{codaSeq.rarefy}},
	\code{\link{codaSeq.propr.aldex.phi}}
}
\examples{
    # data are samples by column, and OTUs by row.
    data(ak_op)
    filt <- codaSeq.filter(ak_op, min.reads=1000, min.prop=0.01, max.prop=1,
        min.occurrence=0.2, samples.by.row=FALSE)
    dim(filt) # should return 179 OTUs and 30 samples

    filt <- codaSeq.filter(ak_op, samples.by.row=FALSE)
    dim(filt) # should return 1061 OTUs and 15 samples
}
