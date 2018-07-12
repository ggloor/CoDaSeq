\name{codaSeq.effect}
\alias{codaSeq.effect}
\title{distributional effect calculation}
\description{
    Returns a vector of the distributional effect size metric, dNEF. This
    is intended for use only when the equivalent metric from ALDEx2 cannot
    be obtained.
}
\usage{
    codaSeq.effect(x, conds, corrected=FALSE)
}
\arguments{
	\item{x}{
		A matrix or dataframe with features by row and samples by columns.
	}
	\item{conds}{
	    A vector with a description of sample labels, of the form
	    conds=c("A","A","B","B"), etc
	}
	\item{corrected}{
	    Adjust to Cohen's d scale or not. Simply multiplies the returned value
	    by the ratio between dNEF and Cohen's d
	}
}
\details{
	Approximates Cohen's d / 1.418 for two normal distributions.
}
\value{
    returns a vector of dNEF
}
\references{
	Please use the citation given by \code{citation(package="CoDaSeq")}
}
\author{
	Greg Gloor, Andrew Fernandes
}
\seealso{
	\code{\link{codaSeq.filter}},
	\code{\link{codaSeq.outlier}},
	\code{\link{codaSeq.rarefy}},
	\code{\link{codaSeq.phi}}
}
\examples{
    # get dataset from ALDEx2 package
    d <- c(rnorm(100,0,1), rnorm(100,2,1))
    e <- c(rnorm(100,2,1), rnorm(100,0,1))
    de <- rbind(d,e)
    conds <- c(rep("A", 100), rep("B",100))

    # convert to clr with an uninformative prior
    # use only the variables with mid-quartile variance as denominator
    # output will have samples by rows
    codaSeq.effect(de, conds, corrected=TRUE)
}
