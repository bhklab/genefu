`compare.proto.cor` <-
function(gene.cor, proto.cor, nn, p.adjust.m=c("none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")) {
	p.adjust.m <- match.arg(p.adjust.m)
	proto <- dimnames(proto.cor)[[1]]
	## select the best two absolute correlations
	best2corix <- t(apply(X=gene.cor, MARGIN=1, FUN=function(x) { return(order(abs(x), decreasing=TRUE)[1:2]) }))
	best2corn <- apply(X=t(apply(X=nn, MARGIN=1, FUN=function(x) { return(x[order(abs(x), decreasing=TRUE)[1:2]]) })), MARGIN=1, FUN=min) ## not perfect since we are not sure that the samples were paired
	best2cor <- cbind(t(apply(X=gene.cor, MARGIN=1, FUN=function(x) { return(x[order(abs(x), decreasing=TRUE)[1:2]]) })), proto.cor[best2corix], best2corn)
	dimnames(best2cor) <- list(dimnames(gene.cor)[[1]], c("r.x1y", "r.x2y", "r.x1x2", "nn"))
	
	rr <- apply(X=best2cor, MARGIN=1, FUN=function(x) { return(cordiff.dep(r.x1y=abs(x[1]), r.x2y=abs(x[2]), r.x1x2=abs(x[3]), n=x[4], alternative="greater")) })
	rr <- rr["p.value",  , drop=TRUE]
	rr <- p.adjust(rr, method=p.adjust.m)
	names(rr) <- dimnames(best2cor)[[1]]
	if(!is.null(names(proto))) {
		rr2 <- data.frame("proto"=names(proto)[best2corix[ ,1]], "cor"=best2cor[ , 1], "signif"=rr, row.names=dimnames(best2cor)[[1]], stringsAsFactors=FALSE)
	} else {
		rr2 <- data.frame("proto"=proto[best2corix[ ,1]], "cor"=best2cor[ , 1], "signif"=rr,  row.names=dimnames(best2cor)[[1]], stringsAsFactors=FALSE)
	}

	return(rr2)
}