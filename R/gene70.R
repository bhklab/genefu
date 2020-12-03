#' @title Function to compute the 70 genes prognosis profile (GENE70) as published by 
#' van't Veer et al. 2002
#'
#' @description
#' This function computes signature scores and risk classifications from gene expression
#'   values following the algorithm used for the 70 genes prognosis profile (GENE70) as 
#'   published by van't Veer et al. 2002.
#'
#' @usage
#' gene70(data, annot, do.mapping = FALSE, mapping,
#'   std = c("none", "scale", "robust"), verbose = FALSE)
#'
#' @param data Matrix of gene expressions with samples in rows and probes in columns, 
#'   dimnames being properly defined.
#' @param annot Matrix of annotations with at least one column named "EntrezGene.ID", 
#'   dimnames being properly defined.
#' @param do.mapping TRUE if the mapping through Entrez Gene ids must be performed (in case 
#'   of ambiguities, the most variant probe is kept for each gene), FALSE otherwise.
#' @param mapping Matrix with columns "EntrezGene.ID" and "probe" used to force the mapping
#'   such that the probes are not selected based on their variance.
#' @param std Standardization of gene expressions: scale for traditional standardization based
#'   on mean and standard deviation, robust for standardization based on the 0.025 and
#'   0.975 quantiles, none to keep gene expressions unchanged.
#' @param verbose TRUE to print informative messages, FALSE otherwise.
#'
#'
#' @return
#' A list with items:
#' - score Continuous signature scores
#' - risk Binary risk classification, 1 being high risk and 0 being low risk.
#' - mapping Mapping used if necessary.
#' - probe If mapping is performed, this matrix contains the correspondence between
#' the gene list (aka signature) and gene expression data
#'
#' @references
#' L. J. van't Veer and H. Dai and M. J. van de Vijver and Y. D. He and A. A. Hart and 
#'   M. Mao and H. L. Peterse and K. van der Kooy and M. J. Marton and A. T. Witteveen and 
#'   G. J. Schreiber and R. M. Kerkhiven and C. Roberts and P. S. Linsley and R. Bernards 
#'   and S. H. Friend (2002) "Gene Expression Profiling Predicts Clinical Outcome of Breast 
#'   Cancer", Nature, 415:530–536.
#'
#' @seealso
#' [genefu::nkis]
#'
#' @examples
#' # load GENE70 signature
#' data(sig.gene70)
#' # load NKI dataset
#' data(nkis)
#' # compute relapse score
#' rs.nkis <- gene70(data=data.nkis)
#' table(rs.nkis$risk)
#' # note that the discrepancies compared to the original publication
#' # are closed to the official cutoff, raising doubts on its exact value.
#' # computation of the signature scores on a different microarray platform
#' # load VDX dataset
#' data(vdxs)
#' # compute relapse score
#' rs.vdxs <- gene70(data=data.vdxs, annot=annot.vdxs, do.mapping=TRUE)
#' table(rs.vdxs$risk)
#'
#' @md
#' @export
if(getRversion() >= "2.15.1")  utils::globalVariables("sig.gene70")

gene70 <-
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
		data <- data[ , intersect(dimnames(sig.gene70)[[1]], dimnames(data)[[2]])]
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

	score <- apply(X=data, MARGIN=1, FUN=function (x, y, method, use) {
	  rr <- NA
    if (sum(complete.cases(x, y)) > 3) {
      rr <- cor(x=x, y=y, method=method, use=use)
    }
    return (rr)
	}, y=sig2[, "average.good.prognosis.profile"], method="spearman", use="complete.obs")
	score <- -score
	official.cutoff <- -0.3
	## cutoff leaving 59% of patients in the poor prognosis group in the original dataset
	risk <- ifelse(score >= official.cutoff, 1, 0)

	names(score) <- names(risk) <- dimnames(data)[[1]]
	
	return(list("score"=score, "risk"=risk, "mapping"=mymapping, "probe"=myprobe))
}
