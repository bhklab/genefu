`ggi` <-
function(data, annot, do.mapping=FALSE, mapping, hg, verbose=FALSE) {

	###########
	#internal functions
	###########
	scale.raw.ggi <- function(ggi, hg) {
	if(length(hg) != length(ggi)) { stop("bad length of hg!") }
			mhg1 <- mean(ggi[hg == 1], na.rm=TRUE) 
			mhg3 <- mean(ggi[hg == 3], na.rm=TRUE)
			mm <- mhg1 + (mhg3 -  mhg1) / 2
			res.scaled <- ((ggi - mm) / (mhg3 - mhg1)) * 2
			return(res.scaled)
	}
	###########
	ggi.gl <- cbind(sig.ggi[ ,c("probe", "EntrezGene.ID")], "coefficient"=ifelse(sig.ggi[ ,"grade"] == 1, -1, 1))
	tt <- sig.score(x=ggi.gl, data=data, annot=annot, do.mapping=do.mapping, mapping=mapping, signed=TRUE, verbose=verbose)
	myprobe <- tt$probe
	mymapping <- tt$mapping
	res <- tt$score
		
	if(!missing(hg)) {
		if(length(hg) != nrow(data)) { stop("hg must have the same length nrow(data)!") }
		mhg1 <- mean(res[hg == 1], na.rm=TRUE) 
		mhg3 <- mean(res[hg == 3], na.rm=TRUE)
		mm <- mhg1 + (mhg3 -  mhg1) / 2
		res.scaled <- ((res - mm) / (mhg3 - mhg1)) * 2
		res <- list("score"=res.scaled, "risk"=ifelse(res.scaled >= 0, 1, 0), "mapping"=mymapping, "probe"=myprobe)	
	} else {
		riskt <- rep(NA, length(res))
		names(riskt) <- names(res)
		res <- list("score"=res, "risk"=riskt, "mapping"=mymapping, "probe"=myprobe)
	}
	
	return (res)
}
