library(scatterplot3d)
library(compositions)
library(zCompositions)
library(ALDEx2)
library(igraph)

source("/Users/ggloor/git/CoDaSeq/chunk/codaSeq_functions.R")
source("/Users/ggloor/git/propr/R/propr-functions.R")
######

d <- read.csv("dols_2016.csv", header=T, row.names=1)

d.tax <- data.frame(d[,1:18])
rownames(d.tax) <- rownames(d)

d.split <- strsplit( as.character(d.tax$species),x)
sp <- vector()

for( i in 1:nrow(d.tax)){
	sp[i] <- paste(d.split[[i]][1],d.split[[i]][2], sep="_")
}



B <- c("TCMID118","TCMID129","TCMID125","TCMID120","TCMID124","TCMID127","TCMID117","TCMID122","TCMID130","TCMID132","TCMID131","TCMID133","TCMID134","TCMID136","TCMID135","TCMID123","TCMID128","TCMID126","TCMID119","TCMID121")
H <- c("TCMID109","TCMID102","TCMID112","TCMID111","TCMID098","TCMID099","TCMID100","TCMID115","TCMID104","TCMID107","TCMID097","TCMID110","TCMID105","TCMID106","TCMID101","TCMID113","TCMID108","TCMID114")

#d.count <- d[,-c(1:18)]
#d.count$TCMID116 <- NULL

d.count <- data.frame(d[,B], d[,H])

d.species <- aggregate(d.count, by=list(sp), FUN=sum)
rownames(d.species) <- d.species$Group.1
d.species$Group.1 <- NULL

d.filt.sparse <- codaSeq.filter(d.species, min.reads=0, min.prop=0,
    min.occurrence=0, samples.by.row=FALSE)

d.n0 <- cmultRepl(t(d.filt.sparse), method="CZM", label=0)

d.clr <- codaSeq.clr(d.n0, samples.by.row=TRUE)

pcx <- prcomp(d.clr)
mv <- sum(pcx$sdev^2)

label.col <- c(rep("red", length(B)), rep("blue", length(H)) )

par(mfrow=c(1,1), mar=c(5,4,4,1))
coloredBiplot(pcx, cex=c(0.5,0.3), varxes=FALSE,
    xlab=paste("PC1: ", round(sum(pcx$sdev[1]^2)/mv, 3)),
    ylab=paste("PC2: ", round(sum(pcx$sdev[2]^2)/mv, 3)),
    scale=0, col="black", xlabs.col=label.col, var.axes=F,
    main=">1% Abundance", cex.lab=1.8, cex.main=1.5)

####
# BARPLOT
####
d.P <- apply(d.n0, 1, function(x){x/sum(x)})
dist.d.clr <- dist(d.clr, method="euclidian")
clust.d <- hclust(dist.d.clr, method="ward.D")

colours <- c("steelblue3", "skyblue1", "indianred1", "mediumpurple1","darkblue", "deepskyblue4","olivedrab3", "#FFED6F", "blue","mediumorchid3", "ivory2", "tan1", "aquamarine3", "#C0C0C0", "royalblue4", "mediumvioletred", "#999933", "#666699", "#CC9933", "#006666", "#3399FF", "#993300", "#CCCC99", "#666666", "#FFCC66", "#6699CC", "#663366", "#9999CC", "#CCCCCC", "#669999", "#CCCC66", "#CC6600", "#9999FF", "#0066CC", "#99CCCC", "#999999", "#FFCC00", "#009999", "#FF9900", "#999966", "#66CCCC", "#339966", "#CCCC33", "#EDEDED")

par(mfrow=c(2,1), mar=c(0,3,2,1)+0.1)
plot(clust.d, main=NULL)
#rect.hclust(clust.d1, 6)
barplot(d.P[,clust.d$order], space=0, col=colours, las=2)

####
# ALDEx
####
conds <- c(rep("B", length(B)), rep("H", length(H)) )
x <- aldex.clr(d.filt.sparse, mc.samples=128)
x.e <- aldex.effect(x, conds)

sig <- rownames(x.e[which(abs(x.e$effect) > 0.8),])
plot(x.e$diff.win, x.e$diff.btw, pch=19, col="grey")
points(x.e[sig,"diff.win"], x.e[sig,"diff.btw"], pch=19, col="red")
abline(0,1)
abline(0,-1)

x.e[sig,]

phi.sma.df <- propr.aldex.phi(x)

d.lo.phi <- subset(phi.sma.df, phi<0.5)

# generate a graphical object
g <- graph.data.frame(d.lo.phi, directed=FALSE)
# get the clusters from the graph object
g.clust <- clusters(g)

# data frame containing the names and group memberships of each cluster
g.df <- data.frame(ID=V(g)$name, cluster=g.clust$membership,
    size=g.clust$csize[g.clust$membership])

g.df$genus <- d.tax[V(g)$name,17]

col=rainbow(max(g.df$cluster))

g.df$col <- col[as.numeric(g.df$cluster)]

g.ids <- V(g)$name
g.names <- d.tax[V(g)$name,17]
V(g)$name <- as.vector(g.names)
V(g)$color <- g.df$col


plot(g, layout=layout.fruchterman.reingold.grid(g, weight=0.05/E(g)$phi),
    vertex.label.color="black", vertex.size=5, vertex.color=V(g)$color)

