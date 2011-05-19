`tamr13` <-
function(data, annot, do.mapping=FALSE, mapping, verbose=FALSE) {

	rest <- resc <- NULL
	for(i in 1:length(sig.tamr13)) {
		rest <- cbind(rest, sig.score(x=sig.tamr13[[i]], data=data, annot=annot, do.mapping=do.mapping, mapping=mapping, signed=TRUE, verbose=verbose)$score)
		resc <- c(resc, attributes(sig.tamr13[[i]])$coefficient)
	}
	dimnames(rest) <- list(dimnames(data)[[1]], names(sig.tamr13))
	names(resc) <- names(sig.tamr13)
	res <- drop(rest %*% resc)
	riskt <- rep(NA, length(res))
	names(riskt) <- names(res)
	res <- list("score"=res, "risk"=riskt)
	
	return (res)
}
