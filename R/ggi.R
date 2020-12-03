#' @title Function to compute the raw and scaled Gene expression Grade Index (GGI)
#'
#' @description
#' This function computes signature scores and risk classifications from gene expression
#'   values following the algorithm used for the Gene expression Grade Index (GGI).
#'
#' @usage
#' ggi(data, annot, do.mapping = FALSE, mapping, hg, verbose = FALSE)
#'
#' @param data Matrix of gene expressions with samples in rows and probes in columns, 
#'   dimnames being properly defined.
#' @param annot	Matrix of annotations with at least one column named "EntrezGene.ID", 
#'   dimnames being properly defined.
#' @param do.mapping TRUE if the mapping through Entrez Gene ids must be performed 
#'   (in case of ambiguities, the most variant probe is kept for each gene), FALSE otherwise.
#' @param mapping Matrix with columns "EntrezGene.ID" and "probe" used to force the 
#'   mapping such that the probes are not selected based on their variance.
#' @param hg Vector containing the histological grade (HG) status of breast cancer 
#'   patients in the dataset.
#' @param verbose TRUE to print informative messages, FALSE otherwise.
#'
#' @return
#' A list with items:
#' - score: Continuous signature scores
#' - risk: Binary risk classification, 1 being high risk and 0 being low risk.
#' - mapping: Mapping used if necessary.
#' - probe: If mapping is performed, this matrix contains the correspondence between 
#' the gene list (aka signature) and gene expression data.
#'
#' @references
#' Sotiriou C, Wirapati P, Loi S, Harris A, Bergh J, Smeds J, Farmer P, Praz V, 
#'   Haibe-Kains B, Lallemand F, Buyse M, Piccart MJ and Delorenzi M (2006) 
#'   "Gene expression profiling in breast cancer: Understanding the molecular basis 
#'   of histologic grade to improve prognosis", Journal of National Cancer Institute, 
#'   98:262–272
#'
#' @seealso
#' [genefu::gene76]
#'
#' @examples
#' # load GGI signature
#' data(sig.ggi)
#' # load NKI dataset
#' data(nkis)
#' # compute relapse score
#' ggi.nkis <- ggi(data=data.nkis, annot=annot.nkis, do.mapping=TRUE,
#'   hg=demo.nkis[ ,"grade"])
#' table(ggi.nkis$risk)
#'
#' @md
#' @export
if(getRversion() >= "2.15.1")  utils::globalVariables("sig.ggi")

ggi <-
function(data, annot, do.mapping=FALSE, mapping, hg, verbose=FALSE) {

	###########
	#internal functions
	###########
	scale.raw.ggi <- function(ggi, hg) {
	if(length(hg) != length(ggi)) { stop("bad length of hg!") }
			mhg1 <- mean(ggi[hg == 1], na.rm=TRUE) 
			mhg3 <- mean(ggi[hg == 3], na.rm=TRUE)
			mm <- mhg1 + (mhg3 -  mhg1) / 2
			res.scaled <- ((ggi - mm) / (mhg3 - mhg1)) * 2
			return(res.scaled)
	}
	###########
	ggi.gl <- cbind(sig.ggi[ ,c("probe", "EntrezGene.ID")], "coefficient"=ifelse(sig.ggi[ ,"grade"] == 1, -1, 1))
	tt <- sig.score(x=ggi.gl, data=data, annot=annot, do.mapping=do.mapping, mapping=mapping, signed=TRUE, verbose=verbose)
	myprobe <- tt$probe
	mymapping <- tt$mapping
	res <- tt$score
		
	if(!missing(hg)) {
		if(length(hg) != nrow(data)) { stop("hg must have the same length nrow(data)!") }
		mhg1 <- mean(res[hg == 1], na.rm=TRUE) 
		mhg3 <- mean(res[hg == 3], na.rm=TRUE)
		mm <- mhg1 + (mhg3 -  mhg1) / 2
		res.scaled <- ((res - mm) / (mhg3 - mhg1)) * 2
		res <- list("score"=res.scaled, "risk"=ifelse(res.scaled >= 0, 1, 0), "mapping"=mymapping, "probe"=myprobe)	
	} else {
		riskt <- rep(NA, length(res))
		names(riskt) <- names(res)
		res <- list("score"=res, "risk"=riskt, "mapping"=mymapping, "probe"=myprobe)
	}
	
	return (res)
}
