\name{codaSeq.PCAplot}
\alias{codaSeq.PCAplot}
\title{PCAplots coloured per group}
\description{
    Generates a PCAplot where the samples can be coloured by group and ellipses can
    be added to show that encompass 75% of the samples.
}
\usage{
	codaSeq.PCAplot(pcx, plot.groups=FALSE, plot.circles=FALSE, plot.loadings=TRUE,
  grp.col=NULL, grp=NULL, main="")
}
\arguments{
	\item{pcx}{
		The output from the prcomp function.
	}
	\item{plot.groups}{
		Plot by group, or not.
	}
	\item{plot.circles}{
		Draw enclosing 75% confidence ellipses or not.
	}
	\item{plot.loadings}{
		Plot loadings (rotations) or not in separate plot.
	}
	\item{grp.col}{
		A vector of colours. One per group.
	}
	\item{grp}{
		A list of vectors with the offsets for the group members.
	}
	\item{main}{
		Plot title
	}
}
\details{
	Provided as an alternative is the coloured biplot function from the compositions R
	package. This is provided as an example approach to add ellipses, etc. to the prcomp
	output plots. It does not perform proper scaling for PCA objects.
	For customization, the user can easily modify the code.
}
\value{
    returns one or two plots showing the loadings and groups in the data.
}
\references{
	Please use the citation given by \code{citation(package="CoDaSeq")}
}
\author{
	Greg Gloor, Jean Macklaim, Andrew Fernandes, Wallace Chan
}
\seealso{
	\code{\link{codaSeq.filter}},
	\code{\link{codaSeq.outlier}},
	\code{\link{codaSeq.rarefy}},
	\code{\link{codaSeq.phi}}
}
\examples{
    # non yet
}
