`sig.score` <-
function(x, data, annot, do.mapping=FALSE, mapping, size=0, cutoff=NA, signed=TRUE, verbose=FALSE) {
	
	if(missing(data) || missing(annot)) { stop("data and annot parameters must be specified") }
	x <- as.data.frame(x, stringsAsFactors=FALSE)
	if(nrow(x) == 0) { stop("empty gene list!"); }

	myprobe <- as.character(x[ ,"probe"])
	mygid <- as.character(x[ ,"EntrezGene.ID"])
	mycoef <- as.numeric(x[ ,"coefficient"])
	names(mycoef) <- names(mygid) <- names(myprobe) <- myprobe

	nix <- order(abs(mycoef), decreasing=TRUE, na.last=NA)
	myprobe <- myprobe[nix]
	mygid <- mygid[nix]
	mycoef <- mycoef[nix]
   
   if(do.mapping) { ## mapping is requested
		gid1 <- mygid
		gid2 <- as.character(annot[ ,"EntrezGene.ID"])
		names(gid2) <- dimnames(annot)[[1]]
		## remove missing and duplicated geneids from the gene list
		rm.ix <- is.na(gid1) | duplicated(gid1)
		gid1 <- gid1[!rm.ix]
	
		rr <- geneid.map(geneid1=gid2, data1=data, geneid2=gid1, verbose=FALSE)
		if(is.na(rr$geneid1[1])) {
			#no gene ids in common
			res <- rep(NA, nrow(data))
			names(res) <- dimnames(data)[[1]]
			gf <- c("mapped"=0, "total"=nrow(x))
			if(verbose) { message(sprintf("probe candidates: 0/%i", nrow(x))) }
			return(list("score"=res, "mapping"=gf, "probe"=cbind("probe"=NA, "EntrezGene.ID"=NA, "new.probe"=NA)))
		}
		nix <- match(rr$geneid2, mygid)
		myprobe <- myprobe[nix]
		mygid <- mygid[nix]
		mycoef <- mycoef[nix]
		gid1 <- rr$geneid2
		if(is.null(names(gid1))) { stop("problem with annotations!") }
		gid2 <- rr$geneid1
		if(is.null(names(gid2))) { stop("problem with annotations!") }
		data <- rr$data1
	
		#change the names of probes in x and data
		names(mycoef) <- names(mygid) <- mygid <- names(myprobe) <- myprobe <- as.character(gid1)
		dimnames(data)[[2]] <- as.character(gid2)
	} else { ## no mapping
		nix <- is.element(myprobe, dimnames(data)[[2]])
		myprobe <- myprobe[nix]
		mygid <- mygid[nix]
		mycoef <- mycoef[nix]
		gid1 <- gid2 <- mygid
		data <- data[ ,myprobe,drop=FALSE]
	}
	if(length(myprobe) == 0) {
		if(verbose) { message(sprintf("probe candidates: 0/%i", size)) }
		tt <- rep(NA, nrow(data))
		names(tt) <- dimnames(data)[[1]]
		return(list("score"=tt, "mapping"=c("mapped"=0, "total"=nrow(x)), "probe"=cbind("probe"=names(gid1), "EntrezGene.ID"=gid1, "new.probe"=names(gid2))))
	}
	
	if(size == 0 || size > nrow(x)) { size <- length(myprobe) }
	nix <- 1:size
	myprobe <- myprobe[nix]
	mygid <- mygid[nix]
	mycoef <- mycoef[nix]
	gid1 <- gid1[nix]
	gid2 <- gid2[nix]
	if(!is.na(cutoff)) {
		nix <- abs(mycoef) > cutoff
		myprobe <- myprobe[nix]
		mygid <- mygid[nix]
		mycoef <- mycoef[nix]
		gid1 <- gid1[nix]
		gid2 <- gid2[nix]
   }
	probe.candp <- myprobe[mycoef >= 0]
	probe.candn <- myprobe[mycoef < 0]
	gf <- length(myprobe)

	gf <- c("mapped"=gf, "total"=nrow(x))
	if(verbose) { message(sprintf("probe candidates: %i/%i",gf[1], gf[2])) }

	nprobe <- c(probe.candp, probe.candn)
	myw <- c("p"=length(probe.candp) / length(nprobe), "n"=length(probe.candn) / length(nprobe))
	res <- rep(0, nrow(data))
	
	if(signed) {
		## consider only the sign of the coefficients
		if(length(probe.candp) > 0) { res <- myw["p"] * (apply(X=data[ ,probe.candp,drop=FALSE], MARGIN=1, FUN=sum, na.rm=TRUE) / apply(X=data[ ,probe.candp,drop=FALSE], MARGIN=1, FUN=function(x) { return(sum(!is.na(x))) })) }
		if(length(probe.candn) > 0) { res <- res - myw["n"] * (apply(X=data[ ,probe.candn,drop=FALSE], MARGIN=1, FUN=sum, na.rm=TRUE) / apply(X=data[ ,probe.candn,drop=FALSE], MARGIN=1, FUN=function(x) { return(sum(!is.na(x))) })) }
	} else {
		## consider the exact value of the coefficients
		if(length(probe.candp) > 0) { res <- myw["p"] * (apply(X=data[ ,probe.candp,drop=FALSE], MARGIN=1, FUN=function(x, y) { nix <- is.na(x); return(sum(x * y, na.rm=TRUE) / sum(y[!nix])) }, y=abs(mycoef[probe.candp]))) }
		if(length(probe.candn) > 0) { res <- res - myw["n"] * (apply(X=data[ ,probe.candn,drop=FALSE], MARGIN=1, FUN=function(x, y) { nix <- is.na(x); return(sum(x * y, na.rm=TRUE) / sum(y[!nix])) }, y=abs(mycoef[probe.candn]))) }
	}
	return(list("score"=res, "mapping"=gf, "probe"=cbind("probe"=names(gid1), "EntrezGene.ID"=gid1, "new.probe"=names(gid2))))
}