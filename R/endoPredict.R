`endoPredict` <-
function(data, annot, do.mapping=FALSE, mapping, verbose=FALSE) {

	## the reference genes are not taken into account due to their absence from most platforms
  #sig2 <- sig.endoPredict
  sig2 <- sig.endoPredict[sig.endoPredict[ , "group"] != "REFERENCE", , drop=FALSE]
	rownames(sig2) <- sig2[ , "probe.affy"]
	gt <- nrow(sig2)
	if(do.mapping) { ## not an affy HGU platform
		gid1 <- as.numeric(as.character(sig2[ ,"EntrezGene.ID"]))
		names(gid1) <- dimnames(sig2)[[1]]
		gid2 <- as.numeric(as.character(annot[ ,"EntrezGene.ID"]))
		names(gid2) <- dimnames(annot)[[1]]
		## remove missing and duplicated geneids from the gene list
		rm.ix <- is.na(gid1) | duplicated(gid1)
		gid1 <- gid1[!rm.ix]
		## mqpping
		rr <- geneid.map(geneid1=gid2, data1=data, geneid2=gid1, verbose=FALSE)
		gm <- length(rr$geneid2)
		mymapping <- c("mapped"=gm, "total"=gt)
		if(!all(is.element(sig2[sig2[ , "group"] == "GOI", "EntrezGene.ID"], rr$geneid1))) { ## if genes of interest are missing
			res <- rep(NA, nrow(data))
			names(res) <- dimnames(data)[[1]]
			if(verbose) { message(sprintf("probe candidates: %i/%i", gm, gt)) }
			return(list("score"=res, "risk"=res, "mapping"=mymapping, "probe"=NA))
		}
		gid1 <- rr$geneid2
		gid2 <- rr$geneid1
		data <- rr$data1
		myprobe <- cbind("probe"=names(gid1), "EntrezGene.ID"=gid1, "new.probe"=names(gid2))
		## change the names of probes in the data
		colnames(data) <- names(gid2) <- names(gid1)
    sig2 <- sig2[colnames(data), , drop=FALSE]
		gm <- ncol(data)
		mymapping <- c("mapped"=gm, "total"=gt)
	} else {
		myprobe <- NA
    nn <- intersect(dimnames(sig2)[[1]], dimnames(data)[[2]])
		data <- data[ , nn]
    sig2 <- sig2[nn, , drop=FALSE]
		gm <- ncol(data)
		mymapping <- c("mapped"=gm, "total"=gt)
	}
	## rename gene names by the gene symbols
	colnames(data) <- rownames(sig2) <- sig2[ , "symbol"]
	
	## transform expressions so they match approximately the scale of Affymetrix data
	if(do.mapping) {
    data <- apply(data, 2, function(x) {
    xx <- (x - quantile(x, probs=0.025, na.rm=TRUE)) / (quantile(x, probs=0.975, na.rm=TRUE) - quantile(x, probs=0.025, na.rm=TRUE)) 
    return((xx * 8) + 6)
  })
  data[!is.na(data) & data < 1] <- 1
  data[!is.na(data) & data > 15] <- 15
}
  
  data <- (data - apply(data, 1, mean, na.rm=TRUE)) + log2(500)
  ## apply transformation factor and offset
  datat <- t(apply(data, 1, function(x, a, b) {
    return((x - b) / a)
  }, a=sig2[ , "a"], b=sig2[ , "b"]))
  data <- matrix(NA, nrow=nrow(data), ncol=ncol(data), dimnames=dimnames(data))
  data[rownames(datat), colnames(datat)] <- datat

	rs <- rs.unscaled <- rsrisk <- rep(NA, nrow(data))
  rs.unscaled <- drop((sig2[ , "weight"] %*% t(data)) - 2.63)
  rs <- sapply(rs.unscaled, function(x) {
    if(!is.na(x)) {
      x <- 1.5 * x + 18.95
      if(x < 0) {
        x <- 0
      } else {
        if(x > 15) {
          x <- 15
        }
      }
    }
    return(x)
  })
  #rsrisk <- ifelse(rs >= 5, 1, 0)
	names(rs) <- names(rs.unscaled) <- names(rsrisk) <- dimnames(data)[[1]]
	return(list("score"=rs, "risk"=rsrisk, "mapping"=mymapping, "probe"=myprobe))
}
