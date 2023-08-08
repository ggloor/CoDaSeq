#' Expected value of Phi from Dirichlet log-ratio distribution
#'
#' @param aldex.clr  An `aldex.clr` S3 object, output by `aldex.clr` from the `ALDEx2` package.
#'
#' @return Returns expected Phi values for pairwise comparisons between
#'   all features.
#'   
#' @details This function is time- and memory-intensive! Requires the `aldex.clr()` 
#' function from the `ALDEx2` package. We ignore all the other measures that
#' are used for trouble-shooting phi.
#' 
#' @export
#'
#' @author Greg Gloor, Jean Macklaim, Wallace Chan
#' 
#' @seealso [codaSeq.clr()], [codaSeq.filter], [codaSeq.rarefy], [codaSeq.outlier]
#' 
#' @examples
#' # load example HMP data from CoDaSeq package
#' data("ak_op")  # feature table: 4347 OTU x 30 samples (15x ak, 15x op)
#' 
#' # CLR using ALDEx2
#' clr.data <- aldex.clr(ak_op,
#'                       conds = c(rep("AK",15),rep("OP",15)),
#'                       denom = "all")
#' 
#' # calculate phi
#' phi<-codaSeq.phi(clr.data)
codaSeq.phi <- function(aldex.clr){

	# calculate expected value of phi
	# a single high phi value will push the component out of consideration
	# a median is right out for memory considerations

	# get first value
	sym.phi <- propr.phisym(t(sapply(getMonteCarloInstances(aldex.clr),
	    function(y){y[,1]})))

	# sum the rest of the values as we proceed through the DIR MC instances
	for(i in 2:numMCInstances(aldex.clr)){
		#print(i)
		sym.phi <- sym.phi + propr.phisym(t(sapply(getMonteCarloInstances(aldex.clr),
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
		col=factor(colnames(sym.phi)[lt.ind[,"col"]]), stringsAsFactors=FALSE)

	#save the lower triangle as an expected value
	sma.df$phi <- sym.phi[lt] /  numMCInstances(aldex.clr)

	return(sma.df)
}
