\name{propr.phisym}
\alias{propr.phisym}
\title{Symmetric Phi Statistics}
\description{
	Returns a matrix where element (i,j) is the symmetric phi statistic between columns i and j of X.
}
\usage{
propr.phismy <- function (X)
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
	Please use the citation given by \code{citation(package=“CoDaSeq”)}
}
\author{
	Greg Gloor, Jean Macklai, Wallace Chan
}
\seealso{
	\code{\link{codaSeq.clr}},
	\code{\link{codaSeq.filter}},
	\code{\link{codaSeq.outlier}},
	\code{\link{codaSeq.rarefy}},
	\code{\link{codaSeq.propr.aldex.phi}}
}
\examples{
	N <- 10
# Number of observations
# Make a data frame with columns a and b roughly proportional
# and columns c and d roughly proportional
 	X <- data.frame(a=(1:N), b=(1:N) * rnorm(N, 10, 0.1),
                 	c=(N:1), d=(N:1) * rnorm(N, 10, 1.0))
 round(propr.phisym(clr(X)),2)
}
