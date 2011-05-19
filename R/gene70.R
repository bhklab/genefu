`gene70` <-
function(data, annot, do.mapping=FALSE, mapping, std=c("none", "scale", "robust"), verbose=FALSE) {
	
	std <- match.arg(std)
	gt <- nrow(sig.gene70)
	if(do.mapping) {
		gid1 <- as.numeric(as.character(sig.gene70[ ,"EntrezGene.ID"]))
		names(gid1) <- dimnames(sig.gene70)[[1]]
		gid2 <- as.numeric(as.character(annot[ ,"EntrezGene.ID"]))
		names(gid2) <- dimnames(annot)[[1]]
		## remove missing and duplicated geneids from the gene list
		rm.ix <- is.na(gid1) | duplicated(gid1)
		gid1 <- gid1[!rm.ix]
	
		rr <- geneid.map(geneid1=gid2, data1=data, geneid2=gid1, verbose=FALSE)
		gm <- length(rr$geneid2)
		if(is.na(rr$geneid1[1])) {
			gm <- 0
			#no gene ids in common
			res <- rep(NA, nrow(data))
			names(res) <- dimnames(data)[[1]]
			gf <- c("mapped"=0, "total"=gt)
			if(verbose) { message(sprintf("probe candidates: 0/%i", gt)) }
			return(list("score"=res, "risk"=res, "mapping"=gf, "probe"=NA))
		}
		gid1 <- rr$geneid2
		gid2 <- rr$geneid1
		data <- rr$data1
		mymapping <- c("mapped"=gm, "total"=gt)
		myprobe <- cbind("probe"=names(gid1), "EntrezGene.ID"=gid1, "new.probe"=names(gid2))
		sig2 <- sig.gene70[names(gid1), , drop=FALSE]
		## change the names of probes in the data
		dimnames(data)[[2]] <- names(gid2) <- names(gid1)
	} else {
		data <- data[ ,intersect(dimnames(sig.gene70)[[1]], dimnames(data)[[2]])]
		sig2 <- sig.gene70[dimnames(data)[[2]], , drop=FALSE]
		gm <- nrow(sig2)
		mymapping <- c("mapped"=gm, "total"=gt)
		myprobe <- NA
	}

	if(verbose && gm != gt) { message(sprintf("%i/%i probes are used to compute the score", gm, gt)) }
	
	## scaling
	switch(std,
	"scale"={
		data <- scale(data, center=TRUE, scale=TRUE)
		if(verbose) { message("standardization of the gene expressions") }
	}, 
	"robust"={
		data <- apply(data, 2, function(x) { return((rescale(x, q=0.05, na.rm=TRUE) - 0.5) * 2) })
		if(verbose) { message("robust standardization of the gene expressions") }
	}, 
	"none"={ if(verbose) { message("no standardization of the gene expressions") } })

	score <- apply(X=data, MARGIN=1, FUN=cor, y=sig2[, "average.good.prognosis.profile"], method="spearman", use="complete.obs")
	score <- -score
	official.cutoff <- -0.3
	## cutoff leaving 59% of patients in the poor prognosis group in the original dataset
	risk <- ifelse(score >= official.cutoff, 1, 0)

	names(score) <- names(risk) <- dimnames(data)[[1]]
	
	return(list("score"=score, "risk"=risk, "mapping"=mymapping, "probe"=myprobe))
}
