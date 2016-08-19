\name{codaSeq.clr}
\alias{codaSeq.clr}
\title{Center Log-Ratio Function.}
\description{
    Returns a matrix of center log-ratio transformed data with samples by row.
    Equivalent to log(x/gx) for every value where gx is the geometic mean of the vector X.
}
\usage{
    codaSeq.clr <- function(x, samples.by.row=TRUE)
}
\arguments{
	\item{x}{
		A matrix or dataframe with samples by rows or columns.
	}
	\item{samples.by.row}{
		TRUE if samples are by row, FALSE if samples are by column.
	}
}
\details{
	Natural log is used for biplots and other exploratory analyses.
}
\value{
    returns a matrix of clr tranformed values with samples in the rows
    and variables in columns
}
\references{
	Please use the citation given by \code{citation(package="CoDaSeq")}
}
\author{
	Greg Gloor, Jean Macklaim, Wallace Chan
}
\seealso{
	\code{\link{codaSeq.filter}},
	\code{\link{codaSeq.outlier}},
	\code{\link{codaSeq.rarefy}},
	\code{\link{codaSeq.propr.aldex.phi}}
}
\examples{
}
