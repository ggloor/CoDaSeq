rdirichlet <- function (n, alpha)
#--------------------------------------------
{
  if(length(n) > 1) n <- length(n)
  #if(length(n) == 0 || as.integer(n) == 0) return(numeric(0))
  #n <- as.integer(n)
  if(n < 0) stop("value(n) can not be negative in rtriang")

  if(is.vector(alpha)) alpha <- t(alpha)
  l <- dim(alpha)[2]
  x <- matrix(rgamma(l * n, t(alpha)), ncol = l, byrow=TRUE)  # Gere le recycling
  return(x / rowSums(x))
}
set.seed(6)

# generate the data, it is just random uniform
s0 <- c(runif(50, min=1e8, max=1e10))
s0[50- <- 2e9
s1 <- c(s0[1:49], 2e10)
s2 <- s1 * runif(length(s1), min=0.8, max=1.2)
s2[50] <- 2e11
s3 <- s1 * runif(length(s1), min=0.8, max=1.2)
s3[50] <- 2e12
s4 <- c(s3[1:49], 2e13)
s <- cbind(s0, s1,as.vector(s2),as.vector(s3), s4)


# proportions with a little bit of multivariate Poisson
# randomness added in
s.p <- apply(s,2,function(x){rdirichlet(1,(x/(0.0001*sum(x))))})
s.clr <- apply(s.p,2, function(x){log10(x) - mean(log10(x))})

colours <- c("steelblue3", "skyblue1", "indianred1", "mediumpurple1","darkblue", "deepskyblue4","olivedrab3", "#FFED6F", "blue","mediumorchid3", "ivory2", "tan1", "aquamarine3", "#C0C0C0", "royalblue4", "mediumvioletred", "#999933", "#666699", "#CC9933", "#006666", "#3399FF", "#993300", "#CCCC99", "#666666", "#FFCC66", "#6699CC", "#663366", "#9999CC", "#CCCCCC", "#669999", "#CCCC66", "#CC6600", "#9999FF", "#0066CC", "#99CCCC", "#999999", "#FFCC00", "#009999", "#FF9900", "#999966", "#66CCCC", "#339966", "#CCCC33", "#EDEDED")

barplot(s.p, space=0, col=colours,names=c(1,2,3,4,5))
pdf("nonlinear_barplot.pdf")
barplot(s.p, space=0, col=colours,names=c(1,2,3,4,5))
dev.off()
