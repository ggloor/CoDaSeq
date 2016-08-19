\name{codaSeq.filter}
\alias{codaSeq.filter}
\title{Filter compositional dataset for 0 values and abundance.}
\description{
	Returns a reduced able of counts where the samples must contain a minimum number of reads,
	and OTUs must be found with a minimum abundance in all remaining samples.
}
\usage{
codaSeq.filter <- function(x, y=tax.vector, min.reads=5000, min.prop=0.001, max.prop=0.025,
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
		The maximum proportional abundance of a read in any sample. Default=0.025.
	}
	\item{min.fraction}{
		The minimum sample proportion of non-0 reads for each variable.
	}
	\item{sample.by.row}{
		True if rows contain samples, false if rows contain variables.
	}
}
\details{
	Filters min/max.prop first, min.fraction second, min/max.prop third.
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
	Greg Gloor, Jean Macklai, Wallace Chan
}
\seealso{
	\code{\link{codaSeq.clr}},
	\code{\link{codaSeq.outlier}},
	\code{\link{codaSeq.rarefy}},
	\code{\link{codaSeq.propr.aldex.phi}}
}
\examples{
}
