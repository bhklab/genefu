#' @title Function to rescale values based on quantiles
#'
#' @description
#' This function rescales values x based on quantiles specified by the user 
#'   such that x' = (x - q1) / (q2 - q1) where q is the specified quantile, 
#'   q1 = q / 2, q2 = 1 - q/2) and x' are the new rescaled values.
#'
#' @usage
#' rescale(x, na.rm = FALSE, q = 0)
#'
#' @param x	
#' @param na.rm	TRUE if missing values should be removed, FALSE otherwise.
#' @param q	Quantile (must lie in [0,1]).
#'
#' @details
#' In order to rescale gene expressions, q = 0.05 yielded comparable scales in 
#'   numerous breast cancer microarray datasets (data not shown).The rational 
#'   behind this is that, in general, 'extreme cases' (e.g. low and high 
#'   proliferation, high and low expression of ESR1, ...) are often present 
#'   in microarray datasets, making the estimation of 'extreme' quantiles 
#'   quite stable. This is specially true for genes exhibiting some 
#'   multi-modality like ESR1 or ERBB2.
#'   
#' @return
#' A vector of rescaled values with two attributes q1 and q1 containing 
#'   the values of the lower and the upper quantiles respectively.
#'
#' @seealso
#' [base::scale]
#'
#' @examples
#' # load VDX dataset
#' data(vdxs)
#' # load NKI dataset
#' data(nkis)
#' # example of rescaling for ESR1 expression
#' par(mfrow=c(2,2))
#' hist(data.vdxs[ ,"205225_at"], xlab="205225_at", breaks=20,
#'   main="ESR1 in VDX")
#' hist(data.nkis[ ,"NM_000125"], xlab="NM_000125", breaks=20,
#'   main="ESR1 in NKI")
#' hist((rescale(x=data.vdxs[ ,"205225_at"], q=0.05) - 0.5) * 2,
#'   xlab="205225_at", breaks=20, main="ESR1 in VDX\nrescaled")
#' hist((rescale(x=data.nkis[ ,"NM_000125"], q=0.05) - 0.5) * 2,
#'   xlab="NM_000125", breaks=20, main="ESR1 in NKI\nrescaled")
#' 
#' @md
#' @export
rescale <-
function(x, na.rm=FALSE, q=0) {
	if(q == 0) {
		ma <- max(x, na.rm=na.rm)
		mi <- min(x, na.rm=na.rm)
	} else {
		ma <- quantile(x, probs=1-(q/2), na.rm=na.rm)
		mi <- quantile(x, probs=q/2, na.rm=na.rm)
	}
	xx <- (x - mi) / (ma - mi)
	attributes(xx) <- list("names"=names(x), "q1"=mi,"q2"=ma)
	return(xx)
}