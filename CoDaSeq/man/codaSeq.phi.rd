\name{codaSeq.phi}
\alias{codaSeq.phi}
\title{Expected Value of Phi From Dirichlet Log-Ratio Distribution}
\description{
	Returns data frame of the lower-triangle of symmetrical phi metric,
	where the value of phi is the expected value of a number of Dirichlet
	Monte-Carlo replicates of the data. This reduces the problem of 0-count and low-count
	features being identified as false positives. These values are
	highly variable because their values range wildly and so the expected
	value is always large.
}
\usage{
codaSeq.phi(aldex.clr)
}
\arguments{
	\item{aldex.clr}{
		Output from aldex.clr function
	}
}
\details{
	Requires aldex.clr function from ALDEx2 Package.
	Param aldex.clr is an S3 object from the aldex.clr function.
	We ignore all the other measures that are used for trouble-shooting phi.
	The function is time and memory intensive
}
\value{
	returns only the expected phi values
}
\references{
	Please use the citation given by \code{citation(package="CoDaSeq")}
}
\author{
	Greg Gloor, Jean Macklaim, Wallace Chan
}
\seealso{
	\code{\link{codaSeq.clr}},
	\code{\link{codaSeq.filter}},
	\code{\link{codaSeq.rarefy}},
	\code{\link{codaSeq.outlier}},
}
\examples{
}
