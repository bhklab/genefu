`compute.proto.cor.meta` <-
function(datas, proto, method=c("pearson", "spearman")) {
	require(survcomp)
	if(!is.list(datas)) {
		if(!all(is.element(proto, dimnames(datas)[[2]]))) { stop("prototypes are not in the dataset!") }
		datasp <- datas[ , proto, drop=FALSE]
		datas <- datas[ , !is.element(dimnames(datas)[[2]], proto)]
		mycor <- matrix(NA, ncol=length(proto), nrow=ncol(datas), dimnames=list(dimnames(datas)[[2]], proto))
		mycorn <- matrix(0, ncol=length(proto), nrow=ncol(datas), dimnames=list(dimnames(datas)[[2]], proto))
		for(i in 1:length(proto)) {
			mycor[ , i] <- apply(X=datas, MARGIN=2, FUN=function(x, y, m) { xx <- NA; if(sum(complete.cases(x, y)) > 1) xx <- cor(x=x, y=y, method=m, use="complete.obs"); return(xx); }, y=datasp[ , i], m=method)
			mycorn[ , i] <- apply(X=datas, MARGIN=2, FUN=function(x, y) { return(sum(complete.cases(x, y))) }, y=datasp[ , i])
		}
	} else {
		nc <- ncol(datas[[1]])
		ncn <- dimnames(datas[[1]])[[2]]
		datast <- datasp <- NULL
		for(k in 1:length(datas)) {
			if(nc != ncol(datas[[k]]) | !all(dimnames(datas[[k]])[[2]] == ncn)) { stop("all the datasets have not the same variables (columns)") }
			if(!all(is.element(proto, dimnames(datas[[k]])[[2]]))) { stop("some prototypes are not in the dataset!") }
			datasp <- c(datasp, list(datas[[k]][ , proto, drop=FALSE]))
			datast <- c(datast, list(datas[[k]][ , !is.element(dimnames(datas[[k]])[[2]], proto)]))
		}
		names(datasp) <- names(datast) <- names(datas)
		datas <- datast
		rm(datast)
		nc <- ncol(datas[[1]])
		ncn <- dimnames(datas[[1]])[[2]]
		
		mycor <- matrix(NA, nrow=nc, ncol=length(proto), dimnames=list(ncn, proto))
		mycorn <- matrix(0, nrow=nc, ncol=length(proto), dimnames=list(ncn, proto))
		for(i in 1:length(proto)) {
			for(j in 1:nc) {
				mycorz <- mycorz.se <- NULL
				nnt <- 0
				for(k in 1:length(datas)) {
					nn <- sum(complete.cases(datas[[k]][ , j], datasp[[k]][ , i]))
					if(nn > 3) {
						mycorz <- c(mycorz, fisherz(cor(x=datas[[k]][ , j], y=datasp[[k]][ , i], method=method, use="complete.obs"), inv=FALSE))
						mycorz.se <- c(mycorz.se, 1/sqrt(nn - 3))
						nnt <- nnt + nn
					} else {
						mycorz <- c(mycorz, NA)
						mycorz.se <- c(mycorz.se, NA)
					}
				}
				mycor[j, i] <- fisherz(combine.est(x=mycorz,x.se=mycorz.se,na.rm=TRUE)$estimate, inv=TRUE)
				mycorn[j, i] <- nnt
			}
		}	
	}
	return(list("cor"=mycor, "cor.n"=mycorn))
}