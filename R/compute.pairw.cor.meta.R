#' @title Function to compute pairwise correlations in a meta-analytical framework
#'
#' @description
#' This function computes meta-estimate of pairwise correlation coefficients for a
#'   set of genes from a list of gene expression datasets.
#'
#' @usage
#' compute.pairw.cor.meta(datas, method = c("pearson", "spearman"))
#'
#' @param datas List of datasets. Each dataset is a matrix of gene expressions with
#'   samples in rows and probes in columns, dimnames being properly defined. All the
#'   datasets must have the same probes.
#' @param method Estimator for correlation coefficient, can be either pearson or spearman.
#'
#' @return
#' A list with items:
#' - cor  Matrix of meta-estimate of correlation coefficients with probes in rows and
#'   prototypes in columns
#' - cor.n Number of samples used to compute meta-estimate of correlation coefficients.
#'
#' @seealso
#' [genefu::map.datasets], [genefu::compute.proto.cor.meta]
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
#' # compute meta-estimate of pairwise correlation coefficients
#' pairwcor <- compute.pairw.cor.meta(datas=datas.mapped$datas, method="pearson")
#' str(pairwcor)
#'
#' @md
#' @importFrom survcomp combine.est fisherz
#' @export
compute.pairw.cor.meta <-
function(datas, method=c("pearson", "spearman")) {
	method <- match.arg(method)
	if(!is.list(datas)) {
		mycor <- cor(x=datas, method=method, use="pairwise.complete.obs")
	} else {
		nc <- ncol(datas[[1]])
		ncn <- dimnames(datas[[1]])[[2]]
		if(length(datas) > 1) {
		    for(k in 2:length(datas)) {
			    if(nc != ncol(datas[[k]]) | !all(dimnames(datas[[k]])[[2]] == ncn)) { stop("all the datasets have not the same variables (columns)") }
		    }
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
						mycorz <- c(mycorz, survcomp::fisherz(cor(x=datas[[k]][ , i], y=datas[[k]][ , j], method=method, use="complete.obs"), inv=FALSE))
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
