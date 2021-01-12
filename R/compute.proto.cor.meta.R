#' @title Function to compute correlations to prototypes in a
#'   meta-analytical framework
#'
#' @description
#' This function computes meta-estimate of correlation coefficients between a set of genes
#'   and a set of prototypes from a list of gene expression datasets.
#'
#' @usage
#' compute.proto.cor.meta(datas, proto, method = c("pearson", "spearman"))
#'
#' @param datas List of datasets. Each dataset is a matrix of gene expressions with samples
#'   in rows and probes in columns, dimnames being properly defined. All the datasets must have the same probes.
#' @param proto	Names of prototypes (e.g. their EntrezGene ID).
#' @param method Estimator for correlation coefficient, can be either pearson or spearman
#'
#' @return
#' A list with items:
#' -cor Matrix of meta-estimate of correlation coefficients with probes in rows and prototypes in columns.
#' -cor.n Number of samples used to compute meta-estimate of correlation coefficients.
#'
#' @seealso
#' [genefu::map.datasets]
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
#'   sample(x=1:nrow(annot.nkis), size=20)))
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
#' str(probecor)
#'
#' @md
#' @importFrom survcomp combine.est fisherz
#' @export
compute.proto.cor.meta <-
function(datas, proto, method=c("pearson", "spearman")) {
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