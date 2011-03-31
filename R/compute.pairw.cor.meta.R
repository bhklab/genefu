`compute.pairw.cor.meta` <-
function(datas, method=c("pearson", "spearman")) {
	require(survcomp)
	if(!is.list(datas)) {
		mycor <- cor(x=datas, method=method, use="pairwise.complete.obs")
	} else {
		nc <- ncol(datas[[1]])
		ncn <- dimnames(datas[[1]])[[2]]
		for(k in 2:length(datas)) {
			if(nc != ncol(datas[[k]]) | !all(dimnames(datas[[k]])[[2]] == ncn)) { stop("all the datasets have not the same variables (columns)") }
		}
		mycor <- matrix(NA, nrow=nc, ncol=nc, dimnames=list(ncn, ncn))
		mycorn <- matrix(0, nrow=nc, ncol=nc, dimnames=list(ncn, ncn))
		for(i in 1:nc) {
			for(j in 1:i) {
				mycorz <- mycorz.se <- NULL
				nnt <- 0
				for(k in 1:length(datas)) {
					if(sum(complete.cases(datas[[k]][ , c(i, j)])) > 1) {
						nn <- sum(complete.cases(datas[[k]][ , c(i, j)]))
						mycorz <- c(mycorz, fisherz(cor(x=datas[[k]][ , i], y=datas[[k]][ , j], method=method, use="complete.obs"), inv=FALSE))
						mycorz.se <- c(mycorz.se, 1/sqrt(nn - 3))
						nnt <- nnt + nn
					} else {
						mycorz <- c(mycorz, NA)
						mycorz.se <- c(mycorz.se, NA)
					}
				}
				mycor[i, j] <- mycor[j, i] <- fisherz(combine.est(x=mycorz,x.se=mycorz.se,na.rm=TRUE)$estimate, inv=TRUE)
				mycorn[i, j] <- mycorn[j, i] <- nnt
			}
		}	
	}
	return(list("cor"=mycor, "cor.n"=mycorn))
}