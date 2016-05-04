# I have placed this script in the directory
# change the path accordingly
source("~/git/propr/R/propr-functions.R")
library(igraph)

###########
# these are functions that you will need. Copy and paste into the console

sma.df <- function(df){
df.cor <- stats::cor(df, use="pairwise.complete.obs")
df.var <- stats::cov(df, use="pairwise.complete.obs")
df.sd <- sqrt(diag(df.var))
r.rf2 <-
    (outer(diag(df.var), diag(df.var), "-")^2 ) /
    (outer(diag(df.var), diag(df.var), "+")^2 - 4 * df.var^2 )
 diag(r.rf2) <- 0
  res.dof     <- nrow(df) - 2
  F           <- r.rf2/(1 - r.rf2) * res.dof
list(b=sign(df.cor) * outer(df.sd, df.sd, "/"),
p=1 - pf(F, 1, res.dof),
r2=df.cor^2)
}

lt.row.min <- function(X){

  result <- numeric(nrow(X) - 1)

  for(i in 2:nrow(X)){

    result[i-1] <- min(X[i, 1:i-1])

  }

  result

}


propr.phisym <- function (X)
{
 Cov    <- stats::var(X)
 tmp    <- 2 * Cov / outer(diag(Cov), diag(Cov), "+")
 return((1-tmp)/(1+tmp))
}
############ end functions
##

# We start with the log-ratio transformed data

# generate the phi statistic
# this is a measure of the concordance of ratios across samples
# low values are better, a good place to start is 0.05, but this has only been
# examined in one dataset
# this is from Lovell's R package in progress
d.n0.sym.phi <- propr.phisym(d.n0.clr)

# generate a symmetrical regression line of best fit
d.n0.sma <- sma.df(d.n0.clr)

# make an empty dataset of the correct size
lt <- which(col(d.n0.sma$b)<row(d.n0.sma$b), arr.ind=FALSE)
lt.ind <- which(col(d.n0.sma$b)<row(d.n0.sma$b), arr.ind=TRUE)

# find the minimum value of each row
d.n0.phi.min <- lt.row.min(d.n0.sym.phi)

# dataframe to hold the info,
# data is a set of vectors where only the lower triangle is kept, because the matrix is symmetrical
# this is needed so the subset function works properly
d.n0.sma.df <- data.frame(row=factor(rownames(d.n0.sma$b)[lt.ind[,"row"]]), col=factor(colnames(d.n0.sma$b)[lt.ind[,"col"]]))
d.n0.sma.df$b <- d.n0.sma$b[lt]
d.n0.sma.df$p <- d.n0.sma$p[lt]
d.n0.sma.df$r2 <- d.n0.sma$r2[lt]
d.n0.sma.df$phi <- d.n0.sym.phi[lt]

# find the set of connections with phi less than some value
# phi.cutoff <- 0.1
d.n0.sma.lo.phi <- subset(d.n0.sma.df, phi < phi.cutoff)

# generate a graphical object
g <- graph.data.frame(d.n0.sma.lo.phi, directed=FALSE)

# overview of all the proportional relationships
# this can take a long time!!!
 plot(g, layout=layout.fruchterman.reingold.grid(g, weight=0.05/E(g)$phi), vertex.size=1, vertex.color="black", vertex.label=NA)

# # get the clusters from the graph object
 g.clust <- clusters(g)
#
# # data frame containing the names and group memberships of each cluster
 g.df <- data.frame(Systematic.name=V(g)$name, cluster=g.clust$membership, cluster.size=g.clust$csize[g.clust$membership])

max.clust <- induced.subgraph(g, which(g.clust$membership %in% which(g.clust$csize == max(g.clust$csize))))
max.clust.names <- V(max.clust)$name

#
# # to see the sized of clusters in the dataset (i.e., the number of members in clusters)
# g.clust$csize
#
# # make a subgraph object for later plotting with a specific cluster size
#
# g.10 <- induced.subgraph(g, which(g.clust$membership %in% which(g.clust$csize == 29)))
# # get the names
# g.10.names <- V(g.10)$name
#
# plot(g.10, layout=layout.fruchterman.reingold.grid(g.10, weight=0.05/E(g.10)$phi), vertex.size=3, vertex.color="white")
#
# # generate a variance variance plot with the clusters overlaid
# # this is the aldex_all.txt data.
# plot(x.all$diff.win, x.all$diff.btw, pch=19,cex=0.3)
# points(x.all[g.10.names,"diff.win"], x.all[g.10.names,"diff.btw"], pch=19, col="red", cex=0.4)
# abline(0,2)
# abline(0,-2)
#
