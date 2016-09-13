\name{ak_op}
\docType{data}
\alias{ak_op}
\title{Selection-based differential sequence variant abundance dataset}
\description{
  This data set is the count table for 15 random op and ak samples from
  the original Human Microbiome dataset. They have been filtered to remove
  any OTU where the read count in all samples is 0.
}
\usage{ak_op}
\format{A dataframe of 3743 OTUs and 30 samples. The first 7 samples are non-selected,
 the last 7 are selected.}
\source{HMP consortium}
\references{
 NIH HMP Working Group et al (2009) Genome Res dx.doi.org/10.1101/gr.096651.109
}
\keyword{datasets}
