# ---- setup ----

# we need these libraries for the colored biplot, and for the 0 replacement function
library(compositions)
library(zCompositions)
source("chunk/codaSeq_functions.R") # commonly used functions

# read the table, with column and row names, columns tab delimited
# samples are by column, variables are by row
d.bf.1 <- read.table("tongue_vs_cheek.txt", header=T, row.names=1, sep="\t")

# move taxonomy info to a vector
tax.0 <- data.frame(d.bf.1$tax)
rownames(tax.0) <- rownames(d.bf.1)

# remove the taxonomy column
d.bf.1$tax <- NULL

# keep only those samples with > min.reads
min.reads <- 5000
# keep only OTUs with an abundance of at least 0.001
min.prop = 0.01
# keep OTUs that are found in at least 30% of samples
cutoff = .3

# using function in codaMB_functions.R
d.subset <- codaSeq.filter(d.bf.1, min.reads=min.reads, min.fraction=cutoff, min.prop=min.prop,
    samples.by.row=FALSE)

tax.subset <- as.data.frame(tax.0[rownames(tax.0) %in% rownames(d.subset),])
rownames(tax.subset) <- rownames(tax.0)[rownames(tax.0) %in% rownames(d.subset)]

# basically I am generating a shortened name for each OTU for display purposes
# you could do this by hand if you don't know grep
# two ways to substitute row and column labels
# gsub replaces any number of characters denoted by the "." that are followed by
# two underscores with nothing from the list of taxa names in the reduced dataset
short_tax <- gsub(".+__", "", tax.subset[,1])

# com is being assigned a vector of values where the "T" is replicated the number of times
# that we observe td_ in the column names of the reduced data subset
com <- c(rep("T", length(grep("td_.", colnames(d.subset))) ),
    rep("B", length(grep("bm_.", colnames(d.subset))) ))
