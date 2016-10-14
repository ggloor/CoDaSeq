\name{codaSeq.stripchart}
\alias{codaSeq.stripchart}
\title{Stripcharts per group}
\description{
    Generates a stripchart where the features are grouped. Examples would include grouping
    OTUs at particular taxonomic levels, or genes into pathways, etc.
}
\usage{
	codaSeq.stripchart(aldex.out=NULL, group.table=NULL, group.label=NULL,
	p.method="wi.eBH", x.axis="effect", effect.cutoff=1, p.cutoff=0, plot.p=TRUE, cex=0.8,
    main=NULL, mar=c(2,12,4,0.5), do.ylab=TRUE, heir=FALSE, heir.base=NULL)
}
\arguments{
	\item{aldex.out}{
		The output of ALDEx2 clr and t.test functions in one data frame.
	}
	\item{group.table}{
		A matrix or data frame with the grouping information. Rownames for this
		and aldex.out must match. Column names are the desired groups.
	}
	\item{group.label}{
		Column to group by from the group.table.
	}
	\item{p.method}{
		Significance test to plot from aldex.out. One of wi.eBH, we.eBH, wi.ep, we.ep.
	}
	\item{x.axis}{
		X axis value to plot from aldex.out. One of effect, diff.btw.
	}
	\item{effect.cutoff}{
		Minimum effect size to color plotted points. Default 1.
	}
	\item{p.cutoff}{
		Maximum adjusted or raw p value to color plotted points. Default 0.
	}
	\item{plot.p}{
		Plot points below the expected raw or adjusted p value.
	}
	\item{cex}{
		size of points to plot.
	}
	\item{main}{
		Plot title
	}
	\item{mar}{
		Left margin size. Adjust manually to find appropriate parameters
	}
	\item{do.ylab}{
		Whether to plot the y axis labels or not. Setting to FALSE gives sets
		yaxt="n".
	}
	\item{heir}{
		Set to FALSE if there is a one-to-one mapping between levels (e.g. OTUs to taxa).
		Set to TRUE if there is a many to one mapping between levels (e.g. genes to functions).
	}
	\item{heir.base}{
		The name or index of the column that contains the same names as the aldex.out
		row names
	}
}
\details{
	This is provided as an example. For customization, modify the code.
	Check to ensure the two input data frames are in the correct format. See documentation.
}
\value{
    returns a series of strip charts clustered by group.
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
