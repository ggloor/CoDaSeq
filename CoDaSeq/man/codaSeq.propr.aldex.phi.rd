\name{propr.aldex.phi}
\alias{pror.aldex.phi}
\title{Expected Value of Phi From Dirichlet Log-Ratio Distribution}
\description{
	Returns data frame of the lower-triangle of symmetrical phi metric,
	where the value of phi is the expected value of a number of Dirichlet
	Monte-Carlo replicates of the data. This reduces the problem of 0-count and low-count
	features being highly variable because their values range wildly and so the expected
	value is always large.
}
\usage{
propr.aldex.phi <- function(aldex.clr)
}
\arguments{
}
\details{
	Requires aldex.clr function from ALDEx2 Package.
	Param aldex.clr is an S3 object from the aldex.clr function.
	We ignore all the other measures that are used for trouble-shooting phi.
	The sma.df function in particular is very time and memory intensive
}
\value{
	\item{sym.phi}{
		Calculated sum of phi values through all DIR MR instances.
	}
	\item{lt}{
		Indice of correct size.
	}
	\item{lt.int}{
		Indice of correct size.
	}
	\item{sma.df}{
		Dataframe to hold info.
	}
	\item{sma.df$phi}{
		Dataframe to hold the lower triangle because the matrix is symmetrical
	}
}
\references{
	Please use the citation given by \code{citation(package=“CoDaSeq”)}
}
\author{
	Greg Gloor, Jean Macklaim, Wallace Chan
}
\seealso{
	\code{\link{codaSeq.clr}},
	\code{\link{codaSeq.filter}},
	\code{\link{codaSeq.rarefy}},
	\code{\link{codaSeq.outlier}},
	\code{\link{codaSeq.propr.phismy}}
}
\examples{
}
