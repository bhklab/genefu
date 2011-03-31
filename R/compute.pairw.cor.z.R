`compute.pairw.cor.z` <-
function(datas, method=c("pearson")) {
	require(survcomp)
	if(!is.list(datas)) {
		mycorz <- mycorz.se <- matrix(NA, nrow=ncol(datas), ncol=ncol(datas), dimnames=list(dimnames(datas)[[2]], dimnames(datas)[[2]]))
		for(i in 1:ncol(datas)) {
			for(j in 1:i) {
				mycorz[i, j] <- mycorz[j, i] <- fisherz(cor(x=datas[ , i], y=datas[ , j], method=method, use="complete.obs"), inv=FALSE)
				mycorz.se[i, j] <- mycorz.se[j, i] <- 1/sqrt(sum(complete.cases(datas[ , c(i, j)])) - 3)
			}
		}
	} else {
		nc <- ncol(datas[[1]])
		ncn <- dimnames(datas[[1]])[[2]]
		mycorz <- mycorz.se <- array(NA, dim=c(nc, nc, length(datas)), dimnames=list(ncn,  ncn,  names(datas)))
		mycorn <- array(0, dim=c(nc, nc, length(datas)), dimnames=list(ncn,  ncn,  names(datas)))
		for(k in 1:length(datas)) {
			if(nc != ncol(datas[[k]])) { stop("all the datasets have not the same number of columns!") }
			for(i in 1:ncol(datas[[k]])) {
				for(j in 1:i) {
					nn <- sum(complete.cases(datas[[k]][ , c(i, j)]))
					mycorz[i, j, k] <- mycorz[j, i, k] <- fisherz(cor(x=datas[[k]][ , i], y=datas[[k]][ , j], method=method, use="complete.obs"), inv=FALSE)
					mycorz.se[i, j, k] <- mycorz.se[j, i, k] <- 1/sqrt(nn - 3)
					mycorn[i, j, k] <- mycorn[j, i, k] <- nn
				}
			}
		}
	}
	return(list("z"=mycorz, "se"=mycorz.se, "nn"=mycorn))
}