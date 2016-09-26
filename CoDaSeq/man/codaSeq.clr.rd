\name{codaSeq.clr}
\alias{codaSeq.clr}
\title{Center Log-Ratio Function.}
\description{
    Returns a matrix of center log-ratio transformed data with samples by row.
    Equivalent to log(x/gx) for every value where gx is the geometic mean of the vector X.
}
\usage{
    codaSeq.clr <- function(x, IQLR=FALSE, samples.by.row=TRUE)
}
\arguments{
	\item{x}{
		A matrix or dataframe with samples by rows or columns.
	}
	\item{IQLR}{
		geometric mean computed on all features if TRUE,
		or on the set of features with variance between the first
		and third quartile. To be used when the data is not centered.
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
    # get dataset from ALDEx2 package
    data(selex)

    # convert to clr with an uninformative prior
    # use only the variables with mid-quartile variance as denominator
    # output will have samples by rows
    selex.clr <- codaSeq.clr(selex + 0.5, IQLR=TRUE, samples.by.row=FALSE)
}
