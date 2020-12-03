#' @title Function to compute the weighted mean and weighted variance of 'x'
#'
#' @description
#' This function allows for computing the weighted mean and weighted variance 
#'   of a vector of continuous values.
#'
#' @usage
#' weighted.meanvar(x, w, na.rm = FALSE)
#'
#' @param x	an object containing the values whose weighted mean is to be computed.
#' @param w	a numerical vector of weights of the same length as x giving the 
#'   weights to use for elements of x.
#' @param na.rm	TRUE if missing values should be removed, FALSE otherwise.
#'
#' @details
#' If w is missing then all elements of x are given the same weight, otherwise 
#'   the weights coerced to numeric by as.numeric. On the contrary of 
#'   weighted.mean the weights are NOT normalized to sum to one. If the sum 
#'   of the weights is zero or infinite, NAs will be returned.
#'
#' @return
#' A numeric vector of two values that are the weighted mean and weighted 
#'   variance respectively.
#'
#' @references
#' http://en.wikipedia.org/wiki/Weighted_variance#Weighted_sample_variance
#'
#' @seealso
#' [stats::weighted.mean]
#' 
#' @examples
#' set.seed(54321)
#' weighted.meanvar(x=rnorm(100) + 10, w=runif(100))
#'
#' @md
#' @export
## weighted mean and weighted variance
## sources:
## http://en.wikipedia.org/wiki/T_test
## http://www.nicebread.de/blog/files/fc02e1635792cb0f2b3cbd1f7e6c580b-10.php
weighted.meanvar <- 
function(x, w, na.rm=FALSE) {
	if(missing(w)) { w <- rep(1, length(x))}
	ii <- complete.cases(x, w)
	if(!na.rm && sum(!ii) > 0) { stop("missing values are present!") } else { 
		w <- as.numeric(w[ii])
		x <- as.numeric(x[ii])
	} 
	sum.w <- sum(w) 
	sum.w2 <- sum(w^2) 
	mean.w <- sum(x * w) / sum(w) 
	var.w <- (sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2, na.rm=na.rm)
	res <- c(mean.w, var.w)
	names(res) <- c("weighted.mean", "weighted.var")
	return(res)
}

