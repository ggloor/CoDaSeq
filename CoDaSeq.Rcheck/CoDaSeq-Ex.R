pkgname <- "CoDaSeq"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('CoDaSeq')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("codaSeq.propr.phismy")
### * codaSeq.propr.phismy

flush(stderr()); flush(stdout())

### Name: propr.phisym
### Title: Symmetric Phi Statistics
### Aliases: propr.phisym

### ** Examples

	N <- 10
# Number of observations
# Make a data frame with columns a and b roughly proportional
# and columns c and d roughly proportional
 	X <- data.frame(a=(1:N), b=(1:N) * rnorm(N, 10, 0.1),
                 	c=(N:1), d=(N:1) * rnorm(N, 10, 1.0))
 round(propr.phisym(clr(X)),2)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
