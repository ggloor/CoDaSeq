\name{codaSeq.rarefy}
\alias{codaSeq.rarefy}
\title{Subsample dataset without replacement.}
\description{
	Returns a table of counts where samples have been sampled without replacement.
	This is included for compatibility, but in general is a bad idea, since it
	results in a loss of information, and distorts the underlying data somewhat.
}
\usage{
codaSeq.rarefy(x, n=1000, samples.by.row=TRUE)
}
\arguments{
	\item{x}{
		A matrix or dataframe containing a count table.
	}
	\item{n}{
		The desired target number of reads per sample. Default=1000.
	}
		\item{samples.by.row}{
		True if rows contain samples, false if rows contain OTUs.
	}
}
\value{
	Returns a matrix in the same orientation as the original with counts per OTU
	reduced to the common sampling depth. This is a constant sum operation.
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
	\code{\link{codaSeq.phi}}
}
\examples{
}
