`intrinsic.cluster.predict` <-
function(sbt.model, data, annot, do.mapping=FALSE, mapping, do.prediction.strength=FALSE, verbose=FALSE) {

	if(missing(data) || missing(annot) || missing(sbt.model)) { stop("data, annot and sbt.mod parameters must be specified") }
	if (!is.matrix(data)) { data <- as.matrix(data) }

	if(is.list(sbt.model)) {
		## retrieve model
		centroids <- sbt.model$centroids
		annot.centroids <- sbt.model$centroids.map
		method.cor <- sbt.model$method.cor
		method.centroids <- sbt.model$method.centroids
		std <- sbt.model$std
		mq <- sbt.model$rescale.q
		mins <- sbt.model$mins
	} else {
		## read model file
		## retrieve centroids
		annot.centroids <- read.csv(sbt.model, comment.char="#", stringsAsFactors=FALSE)
		dimnames(annot.centroids)[[1]] <- annot.centroids[ , "probe"]
		## retrieve model parameters
		rr <- readLines(con=sbt.model)[-1]
		rr <- rr[sapply(rr, function(x) { return(substr(x=x, start=1, stop=1) == "#") })]
		nn <- unlist(lapply(X=rr, FUN=function(x) { x <- unlist(strsplit(x=unlist(strsplit(x=x, split=":"))[1], split=" ")); x <- x[length(x)]; return(x); }))
		rr2 <- unlist(lapply(X=rr[is.element(nn, c("method.cor", "method.centroids", "std", "rescale.q", "mins"))], FUN=function(x) { x <- unlist(strsplit(x=unlist(strsplit(x=x, split=":"))[2], split=" ")); x <- x[length(x)]; return(x); }))
		method.cor <- rr2[nn == "method.cor"]
		method.centroids <- rr2[nn == "method.centroids"]
		std <- rr2[nn == "std"]
		mq <- as.numeric(rr2[nn == "rescale.q"])
		mins <- as.numeric(rr2[nn == "mins"])
		rr <- rr[!is.element(nn, c("method.cor", "method.centroids", "std", "rescale.q", "mins"))]
		nn <- nn[!is.element(nn, c("method.cor", "method.centroids", "std", "rescale.q", "mins"))]
		cent <- lapply(X=rr, FUN=function(x) { x <- as.numeric(unlist(strsplit(x=unlist(strsplit(x=x, split=":"))[2], split=" "))); x <- x[!is.na(x)]; return(x)})
		centroids <- NULL
		for(i in 1:length(cent)) { centroids <- cbind(centroids, cent[[i]]) }
		nn <- unlist(lapply(X=strsplit(x=nn, split="centroid."), FUN=function(x) { return(x[[2]]) }))
		dimnames(centroids) <- list(dimnames(annot.centroids)[[1]], nn)
	}

	number.cluster <- ncol(centroids)
	if(is.null(dimnames(centroids)[[2]])) { name.cluster <- paste("cluster", 1:ncol(centroids), sep=".") } else { name.cluster <- dimnames(centroids)[[2]] }
	
	gt <- nrow(centroids)
	#mapping
	if(do.mapping) {
		#mapping through EntrezGene.ID
		#remove (arbitrarily) duplicated gene ids
		centroids.gid <- as.character(annot.centroids[ ,"EntrezGene.ID"])
		names(centroids.gid) <- as.character(annot.centroids[ , "probe"])
		myx <- !duplicated(centroids.gid) & !is.na(centroids.gid)
		centroids.gid <- centroids.gid[myx]
		annot.centroids <- annot.centroids[myx, , drop=FALSE]
		centroids <- centroids[myx, , drop=FALSE]
		gid <- as.character(annot[ ,"EntrezGene.ID"])
		names(gid) <- as.character(annot[ ,"probe"])
		if(missing(mapping)) { ## select the most variant probes using annotations
			## if multiple probes for one gene, keep the most variant
			rr <- geneid.map(geneid1=gid, data1=data, geneid2=centroids.gid, verbose=FALSE)
			nn <- match(rr$geneid2, centroids.gid)
			nn <- nn[!is.na(nn)]
			centroids.gid <- centroids.gid[nn]
			annot.centroids <- annot.centroids[nn, ]
			centroids <- centroids[nn, , drop=FALSE]
			data <- rr$data1
		} else { # use a predefined mapping
			nn <- as.character(mapping[ , "EntrezGene.ID"])
			# keep only centroids genes with mapped probes
			myx <- is.element(centroids.gid, nn)
			centroids.gid <- centroids.gid[myx]
			annot.centroids <- annot.centroids[myx, , drop=FALSE]
			centroids <- centroids[myx, , drop=FALSE]
			pp <- as.character(mapping[match(centroids.gid, nn), "probe"])
			myx <- is.element(pp, dimnames(data)[[2]])
			centroids.gid <- centroids.gid[myx]
			annot.centroids <- annot.centroids[myx, , drop=FALSE]
			centroids <- centroids[myx, , drop=FALSE]
			pp <- pp[myx]
			data <- data[ , pp, drop=FALSE]
		}
	}
	else {
		if(all(!is.element(dimnames(data)[[2]], dimnames(centroids)[[1]]))) { stop("no probe in common -> annot or mapping parameters are necessary for the mapping process!") }
		## no mapping are necessary
		myx <- intersect(dimnames(centroids)[[1]], dimnames(data)[[2]])
		data <- data[ ,myx, drop=FALSE]
		centroids <- centroids[myx, , drop=FALSE]
	}
	centroids.map <- cbind("probe"=dimnames(data)[[2]], "probe.centroids"=dimnames(centroids)[[1]], "EntrezGene.ID"=as.character(annot[dimnames(data)[[2]], "EntrezGene.ID"]))
	dimnames(centroids.map)[[1]] <- dimnames(data)[[2]]
	gm <- nrow(centroids)
	if(gm == 0 || (sum(is.na(data)) / length(data)) > 0.9) { ## no mapping or too many missing values
		ncl <- rep(NA, nrow(data))
		names(ncl) <- dimnames(data)[[1]]
		nproba <- ncor <- matrix(NA, nrow=nrow(data), ncol=ncol(centroids), dimnames=list(dimnames(data)[[1]], name.cluster))
		ps.res <- NULL
		if(do.prediction.strength) { ps.res <- list("ps"=NA, "ps.cluster"=ncor[ , 1], "ps.individual"=ncl)}
		tt <- matrix(NA, ncol=nrow(centroids.map), nrow=nrow(data), dimnames=list(dimnames(data)[[1]], dimnames(centroids.map)[[1]]))
		return(list("subtype"=ncl, "subtype.proba"=nproba, "cor"=ncor, "prediction.strength"=ps.res, "centroids.map"=centroids.map, "profiles"=tt))
	}
	if(verbose) { message(sprintf("%i/%i probes are used for clustering", gm, gt)) }

	#standardization of the gene expressions
	switch(std,
	"scale"={
		data <- scale(data, center=TRUE, scale=TRUE)
		if(verbose) { message("standardization of the gene expressions") }
	}, 
	"robust"={
		data <- apply(data, 2, function(x) { return((rescale(x, q=mq, na.rm=TRUE) - 0.5) * 2) })
		if(verbose) { message("robust standardization of the gene expressions") }
	}, 
	"none"={ if(verbose) { message("no standardization of the gene expressions") } })
	
	## apply the nearest centroid classifier to classify the samples again
	ncor <- t(apply(X=data, MARGIN=1, FUN=function(x, y, method.cor) { return(cor(x, y, method=method.cor, use="complete.obs")) }, y=centroids, method.cor=method.cor))
	#nproba <- t(apply(X=ncor, MARGIN=1, FUN=function(x) { return(abs(x) / sum(abs(x), na.rm=TRUE)) }))
	## negative correlations are truncated to zero since they have no meaning for subtypes identification
	nproba <- t(apply(X=ncor, MARGIN=1, FUN=function(x) { x[x < 0] <- 0; return(x / sum(x, na.rm=TRUE)); }))
	dimnames(ncor) <- dimnames(nproba) <- list(dimnames(data)[[1]], name.cluster)
	ncl <- apply(X=ncor, MARGIN=1, FUN=function(x) { return(order(x, decreasing=TRUE, na.last=TRUE)[1]) })
	names(ncl) <- dimnames(data)[[1]]
	## names of identified clusters
	ncln <- name.cluster[ncl]
	names(ncln) <- dimnames(data)[[1]]
	
	## if one or more subtypes have not been identified, remove them for prediction strength
	myx <- sort(unique(ncl))
	myx <- myx[!is.na(myx)]
	name.cluster2 <- name.cluster[myx]
	number.cluster2 <- length(myx)
	
	ps.res <- ncl2 <- NULL
	if(do.prediction.strength) {
		## compute the clustering and cut the dendrogram
		## hierarchical clustering with correlation-based distance and average linkage
		hcl <- amap::hcluster(x=data, method="correlation", link="average")
		mins.ok <- stop.ok <- FALSE
		nbc <- number.cluster2
		nclust.best <- 1
		while(!mins.ok && !stop.ok) { ## until each cluster contains at least mins samples
			cl <- cutree(tree=hcl, k=nbc)
			tt <- table(cl)
			if(sum(tt >= mins) >= number.cluster2) {
				if(nbc > number.cluster2) { ## put NA for clusters with less than mins samples
					td <- names(tt)[tt < mins]
					cl[is.element(cl, td)] <- NA
					## rename the clusters
					ucl <- sort(unique(cl))
					ucl <- ucl[!is.na(ucl)]
					cl2 <- cl
					for(i in 1:number.cluster2) { cl2[cl == ucl[i] & !is.na(cl)] <- i }
					cl <- cl2
				}
				nclust.best <- number.cluster2
				mins.ok <- TRUE
			} else {
				if(sum(tt >= mins) > nclust.best) {
					nbc.best <- nbc
					nclust.best <- sum(tt >= mins)
				}
				nbc <- nbc + 1
				if(nbc > (nrow(data) - (number.cluster2 * mins))) {
					warning(sprintf("impossible to find %i main clusters with at least %i individuals!", number.cluster2, mins))
					stop.ok <- TRUE
				}
			}
			if(stop.ok) { ## no convergence for the clustering with mininmum set of individuals
				cl <- cutree(tree=hcl, k=nbc.best)
				tt <- table(cl)
				td <- names(tt)[tt < mins]
				cl[is.element(cl, td)] <- NA
				## rename the clusters
				ucl <- sort(unique(cl))
				ucl <- ucl[!is.na(ucl)]
				cl2 <- cl
				for(i in 1:nclust.best) { cl2[cl == ucl[i] & !is.na(cl)] <- i }
				cl <- cl2
			}
		}
		## compute the centroids
		## take the core samples in each cluster to compute the centroid
		## not feasible due to low intra correlation within clusters!!!
		## minimal pairwise cor of approx 0.3
		#cl2 <- cutree(tree=hcl, h=0.7)
		#table(cl, cl2) to detect which core cluster of samples for which cluster.
		cl.centroids <- matrix(NA, nrow=ncol(data), ncol=nclust.best, dimnames=list(dimnames(data)[[2]], paste("cluster", 1:nclust.best, sep=".")))
		for(i in 1:ncol(cl.centroids)) {
			switch(method.centroids, 
			"mean"={ cl.centroids[ ,i] <- apply(X=data[cl == i & !is.na(cl), ,drop=FALSE], MARGIN=2, FUN=mean, na.rm=TRUE, trim=0.025) }, 
			"median"={ cl.centroids[ ,i] <- apply(X=data[cl == i & !is.na(cl), ,drop=FALSE], MARGIN=2, FUN=median, na.rm=TRUE) }, 
			"tukey"={ cl.centroids[ ,i] <- apply(X=data[cl == i & !is.na(cl), ,drop=FALSE], MARGIN=2, FUN=tbrm, na.rm=TRUE, C=9) })
		}
		#apply the nearest centroid classifier to classify the samples again
		ncor2 <- t(apply(X=data, MARGIN=1, FUN=function(x, y, z) { return(cor(x, y, method=z, use="complete.obs")) }, y=cl.centroids, z=method.cor))
		nproba2 <- t(apply(X=ncor2, MARGIN=1, FUN=function(x) { return(abs(x) / sum(abs(x), na.rm=TRUE)) }))
		dimnames(ncor2) <- dimnames(nproba2) <- list(dimnames(data)[[1]], dimnames(cl.centroids)[[2]])
		ncl2 <- apply(X=ncor2, MARGIN=1, FUN=function(x) { return(order(x, decreasing=TRUE)[1]) })
		names(ncl2) <- dimnames(data)[[1]]
		## rename clusters since we do not expect to get the same id per cluster
		## this avoids a warning in ps.cluster
		uncl <- sort(unique(ncl))
		uncl <- uncl[!is.na(uncl)]
		nclt <- ncl
		for(mm in 1:length(uncl)) {
			nclt[ncl == uncl[mm]] <- mm
		}
		uncl2 <- sort(unique(ncl2))
		uncl2 <- uncl2[!is.na(uncl2)]
		ncl2t <- ncl2
		for(mm in 1:length(uncl2)) {
			ncl2t[ncl2 == uncl2[mm]] <- mm
		}
		#prediction strength
		ps.res <- ps.cluster(cl.tr=ncl2t, cl.ts=nclt, na.rm=TRUE)
		## put NA for clusters which are potentially not present in the dataset
		tt <- rep(NA, length(name.cluster))
		names(tt) <- name.cluster
		tt[name.cluster2] <- ps.res$ps.cluster
		ps.res$ps.cluster <- tt
	}
	
	
	return(list("subtype"=ncln, "subtype.proba"=nproba, "cor"=ncor, "prediction.strength"=ps.res, "subtype.train"=ncl2, "profiles"=data, "centroids.map"=centroids.map))
}