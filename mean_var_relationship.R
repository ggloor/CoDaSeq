library(zCompositions)
library(ALDEx2)
source("/Users/ggloor/git/CoDaSeq/chunk/codaSeq_functions.R")

meta <- read.table("metadata.txt", header=T, row.names=1, check.names=F)

# Group:sample:lane
nms <- paste(meta[,2], meta[,3], meta[,1], sep=":")

d <- read.table("~/git/9545/R_section/countfinal2.file", header=T, row.names=1, check.names=F)

# double check that the column names of d and rownames of meta
# are congruent - they are

# change the column names to something more informative
colnames(d) <- nms
d.wt1 <-  d[which(apply(d[,1:7], 1, min) > 0),1:7]
# mean, var if use ratios on simplex
d.n0 <- cmultRepl(t(d.wt1), label=0, method="CZM", output="counts")
d.n0 <- t(d.wt1)
clr.d <- apply(d.n0, 1, function(x){log(x) - mean(log(x))})
mean.clr <- apply(clr.d, 1, mean)
var.clr <- apply(clr.d, 1, var)

# mean, var if use absolute numbers on simplex
mean.abs <- apply(d.n0, 2, mean)
var.abs <- apply(d.n0, 2, var)

x <- aldex.clr(d.wt1, mc.samples=1028)

cl2p <- NULL
sum <- vector("list", 2)

for (m in getMonteCarloInstances(x)){
   cl2p <- cbind(cl2p,m)
}

sum$mean <- t(apply(cl2p, 1, mean))
sum$var <- t(apply(cl2p, 1, var))
nms <- getFeatureNames(x)

pdf("mean_var.pdf", height=10, width=4.5)
par(mfrow=c(2,1), mar=c(4,5,3,1))
plot(mean.abs,var.abs , pch=19, cex.main=1.8, cex.lab=1.8, log="xy", col=rgb(0,0,0,0.1), main="Abs Pt")
plot(mean.abs,var.clr , pch=19, cex.main=1.8, cex.lab=1.8, log="x", col=rgb(0,0,0,0.1), main="Ratio Pt")
dev.off()

pdf("mean_var_Bayes.pdf", height=10, width=4.5)
par(mfrow=c(2,1), mar=c(4,5,3,1))
plot(mean.abs,var.clr , pch=19, cex.main=1.8, cex.lab=1.8, log="x", col=rgb(0,0,0,0.1), main="Ratio Pt")
points(mean.abs[nms], sum$var , pch=19, xlab="mean.abs", ylab="var.clr", cex.main=1.8, log="x",cex.lab=1.8, col=rgb(1,0,0,0.01), main="Ratio Bay")
dev.off()
