\name{codaSeq.PCAplot}
\alias{codaSeq.PCAplot}
\title{PCAplots coloured per group}
\description{
    Generates a PCAplot where the samples can be coloured by group and ellipses can
    be added to show that encompass 75% of the samples. The loadings can also be
    colored. Either samples or loadings can be colored, but not both.
}
\usage{
	codaSeq.PCAplot(pcx, plot.groups=FALSE, plot.circles=FALSE, plot.loadings=TRUE,
  grp.col=NULL, grp=NULL, loadings.grp=NULL, loadings.col=NULL, loadings.sym=19, main="")
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
		Plot loadings (rotations) or not in the same plot.
	}
	\item{grp.col}{
		A vector of colours. One per group.
	}
	\item{grp}{
		A list of vectors with the offsets for the group members.
	}
	\item{loadings.grp}{
		A list of loading vectors with the offsets for the group of loadings.
	}
	\item{loadings.col}{
		A vector of loadings colours. One per group.
	}
	\item{loadings.sym}{
	    A value or vector of loadings symbols. Default is 19
	    Any other value must be a vector that is the length of the color vector
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
    returns one plot showing the loadings and groups in the data.
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
