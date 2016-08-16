\name{codaSeq.outlier}
\alias{codaSeq.outlier}
\title{Identifying sample Outliers.}
\description{
	Returns a list of proportional contribution to group variance,
	sample names that are outliers, and sample names that are not outliers.
}
\usage{
codaSeq.outlier <- function(x, plot.me=TRUE, col=rgb(1,0,0,0.3))
}
\arguments{
	\item{x}{
		A matrix or data frame with clr transformed values, with sample by row. 
	}
	\item{plot.me}{
		A logical value determining if a histogram should be plotted of the variance contribution per sample.
	} 
	\item{col}{
		RGB values for your selection of colour.
	}
}
\details{
	Samples must be grouped.
	This approach makes no sense across groups. If you do not know if you have natural groups,
	ignore this step and exam your data by PCA.
	Outliers are defined as those contributing greater than the median plus twice the 
	interquartile range of the sample variance to the total.
}
\value{
	Returns list
	\item{sample.var}{
		Proportional variance contributions for each sample.
	}
	\item{bad}{
		Rownames of outlier samples.
	}
	\item{good}{
		Rownames of non-outlier samples.
	}
}
\references{
	Please use the citation given by \code{citation(package=“CoDaSeq”)}
}
\author{
	Greg Gloor, Jean Macklai, Wallace Chan
}
\seealso{
	\code{\link{codaSeq.clr}},
	\code{\link{codaSeq.filter}},
	\code{\link{codaSeq.rarefy}},
	\code{\link{codaSeq.propr.phismy}},
	\code{\link{codaSeq.propr.aldex.phi}}
}
\examples{
}
	 

	
	