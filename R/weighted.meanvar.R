## weighted mean and weighted variance
## sources:
## http://en.wikipedia.org/wiki/T_test
## http://www.nicebread.de/blog/files/fc02e1635792cb0f2b3cbd1f7e6c580b-10.php
`weighted.meanvar` <- 
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

