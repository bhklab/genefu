#' @title Function to statistically compare correlation to prototypes
#'
#' @description
#' This function performs a statistical comparison of the correlation 
#'  coefficients as computed between each probe and prototype.
#'
#' @usage
#' compareProtoCor(gene.cor, proto.cor, nn,
#'  p.adjust.m = c("none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"))
#'
#' @param gene.cor Correlation coefficients between the probes and each of the prototypes.
#' @param proto.cor Pairwise correlation coefficients of the prototypes.
#' @param nn Number of samples used to compute the correlation coefficients between
#'   the probes and each of the prototypes.
#' @param p.adjust.m Correction method as defined in p.adjust.
#'
#'
#' @return
#' Data frame with probes in rows and with three columns: 
#'  "proto" is the prototype to which the probe is the most correlated,  
#'  "cor" is the actual correlation, and "signif" is the (corrected) p-value
#'  for the superiority of the correlation to this prototype compared to the 
#'  second highest correlation.
#'
#' @seealso
#' [genefu::compute.proto.cor.meta], [genefu::compute.pairw.cor.meta]
#'
#' @examples
#' # load VDX dataset
#' data(vdxs)
#' # load NKI dataset
#' data(nkis)
#' # reduce datasets
#' ginter <- intersect(annot.vdxs[ ,"EntrezGene.ID"], annot.nkis[ ,"EntrezGene.ID"])
#' ginter <- ginter[!is.na(ginter)][1:30]
#' myx <- unique(c(match(ginter, annot.vdxs[ ,"EntrezGene.ID"]),
#'   sample(x=1:nrow(annot.vdxs), size=20)))
#' data2.vdxs <- data.vdxs[ ,myx]
#' annot2.vdxs <- annot.vdxs[myx, ]
#' myx <- unique(c(match(ginter, annot.nkis[ ,"EntrezGene.ID"]),
#' sample(x=1:nrow(annot.nkis), size=20)))
#' data2.nkis <- data.nkis[ ,myx]
#' annot2.nkis <- annot.nkis[myx, ]
#' # mapping of datasets
#' datas <- list("VDX"=data2.vdxs,"NKI"=data2.nkis)
#' annots <- list("VDX"=annot2.vdxs, "NKI"=annot2.nkis)
#' datas.mapped <- map.datasets(datas=datas, annots=annots, do.mapping=TRUE)
#' # define some prototypes
#' protos <- paste("geneid", ginter[1:3], sep=".")
#' # compute meta-estimate of correlation coefficients to the three prototype genes
#' probecor <- compute.proto.cor.meta(datas=datas.mapped$datas, proto=protos,
#'   method="pearson")
#' # compute meta-estimate of pairwise correlation coefficients between prototypes
#' datas.proto <- lapply(X=datas.mapped$datas, FUN=function(x, p) {
#'   return(x[ ,p,drop=FALSE]) }, p=protos)
#' protocor <- compute.pairw.cor.meta(datas=datas.proto, method="pearson")
#' # compare correlation coefficients to each prototype
#' res <- compareProtoCor(gene.cor=probecor$cor, proto.cor=protocor$cor,
#' nn=probecor$cor.n, p.adjust.m="fdr")
#' head(res)
#'
#' @md
#' @export
compareProtoCor <-
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
