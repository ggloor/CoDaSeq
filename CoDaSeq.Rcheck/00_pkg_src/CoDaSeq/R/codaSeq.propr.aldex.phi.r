
codaSeq.propr.aldex.phi <- function(aldex.clr){

	# calculate expected value of phi
	# a single high phi value will push the component out of consideration
	# a median is right out for memory considerations

	# get first value
	sym.phi <- codaSeq.propr.phisym(t(sapply(getMonteCarloInstances(aldex.clr),
	    function(y){y[,1]})))

	# sum the rest of the values as we proceed through the DIR MC instances
	for(i in 2:numMCInstances(aldex.clr)){
		#print(i)
		sym.phi <- sym.phi + codaSeq.propr.phisym(t(sapply(getMonteCarloInstances(aldex.clr),
		    function(y){y[,i]})))
	}
	##### Done ALDEx2 stuff

	# make indices of the correct size
	lt <- which(col(sym.phi)<row(sym.phi), arr.ind=FALSE)
	lt.ind <- which(col(sym.phi)<row(sym.phi), arr.ind=TRUE)

	# dataframe to hold the info,
	# data is a set of vectors where only the lower triangle is kept, because the matrix
	#    is symmetrical
	# this is needed so subsequent subset function works properly
	sma.df <- data.frame(row=factor(rownames(sym.phi)[lt.ind[,"row"]]),
		col=factor(colnames(sym.phi)[lt.ind[,"col"]]))

	#save the lower triangle as an expected value
	sma.df$phi <- sym.phi[lt] /  numMCInstances(aldex.clr)

	return(sma.df)
}
