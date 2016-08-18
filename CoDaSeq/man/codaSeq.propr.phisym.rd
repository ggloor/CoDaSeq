\name{codaSeq.propr.phisym}
\alias{codaSeq.propr.phisym}
\title{Symmetric Phi Statistics}
\description{
	Returns a matrix where element (i,j) is the symmetric phi statistic between columns i and j of X.
}
\usage{
codaSeq.propr.phismy <- function (X)
}
\arguments{
	\item{X}{
		A matrix or data frame of centered log ratio transformation.
	}
}
\details{
	X should be the result of a centered log-ratio transformation.
}
\value{
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
	\code{\link{codaSeq.outlier}},
	\code{\link{codaSeq.rarefy}},
	\code{\link{codaSeq.propr.aldex.phi}}
}
\examples{
}
