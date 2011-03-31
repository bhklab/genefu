`bimod` <-
function(x, data, annot, do.mapping=FALSE, mapping, model=c("E", "V"), do.scale=TRUE, verbose=FALSE, ...) {
	require(mclust)
	model <- match.arg(model)
	dd <- sig.score(x=x, data=data, annot=annot, do.mapping=do.mapping, mapping=mapping, verbose=verbose, ...)$score
	if(do.scale) { dd <- (rescale(x=dd, q=0.05, na.rm=TRUE) - 0.5) * 2 }
	cc.ix <- complete.cases(dd)

	mystatus <- mystatus.proba <- rep(NA, nrow(data))
	names(mystatus) <- names(mystatus.proba) <- dimnames(data)[[1]]
	res <- matrix(NA, nrow=3, ncol=2, dimnames=list(c("mean", "variance", "proportion"), paste("cluster", 1:2, sep=".")))
	mybic <- matrix(NA, nrow=10, ncol=1, dimnames=list(1:10, model))
	
	if(sum(cc.ix) >= 10) {	
		#How many Gaussians?
		rr <- Mclust(data=dd[cc.ix], modelNames=model, G=1:10)
		oo <- order(rr$BIC, decreasing=TRUE)[1]
		if(oo != 2) { warning(sprintf("%i is the most likely number of Gaussians!", oo)) }
		mybic <- rr$BIC

		#Only 2 Gaussians
		rr2 <- Mclust(data=dd[cc.ix], modelNames=model, G=2)
		if(is.null(rr2[[1]])) { ## EM algorithm did not converge
			return(list("status"=mystatus, "status1.proba"=mystatus.proba, "gaussians"=res, "BIC"=rr$BIC, "x"=dd))
		}
		res[1, ] <- rr2$parameters$mean
		res[2, ] <- rr2$parameters$variance$sigmasq
		res[3, ] <- rr2$parameters$pro
		
		## bimodality index (BI)
		smd <- abs(res[1, 2] - res[1, 1]) / sqrt((res[2, 2]^2 + res[2, 1]^2) / 2)
		bi <- sqrt(res[3, 2] * (1 - res[3, 2])) * smd

		#classification
		mystatus[cc.ix] <- as.numeric(rr2$classification == 2)
		mystatus.proba[cc.ix] <- rr2$z[ , 2, drop=TRUE]
	}
	return(list("status"=mystatus, "status1.proba"=mystatus.proba, "gaussians"=res, "BIC"=mybic,  "BI"=bi, "x"=dd))
}