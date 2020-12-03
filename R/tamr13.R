#' @title Function to compute the risk scores of the tamoxifen resistance 
#' signature (TAMR13)
#'
#' @description
#' This function computes signature scores from gene expression values 
#'   following the algorithm used for the Tamoxifen Resistance signature (TAMR13).
#'
#' @usage
#' tamr13(data, annot, do.mapping = FALSE, mapping, verbose = FALSE)
#'
#' @param data Matrix of gene expressions with samples in rows and probes 
#'   in columns, dimnames being properly defined.
#' @param annot Matrix of annotations with at least one column named 
#'   "EntrezGene.ID", dimnames being properly defined.
#' @param do.mapping TRUE if the mapping through Entrez Gene ids must be 
#'   performed (in case of ambiguities, the most variant probe is kept for 
#'   each gene), FALSE otherwise.
#' @param mapping Matrix with columns "EntrezGene.ID" and "probe" used to 
#' force the mapping such that the probes are not selected based on their variance.
#' @param verbose TRUE to print informative messages, FALSE otherwise.
#'
#' @return
#' A list with items:
#' - score: Continuous signature scores.
#' - risk: Binary risk classification, 1 being high risk and 0 being low 
#'   risk (not implemented, the function will return NA values).
#'
#' @references
#' Loi S, Haibe-Kains B, Desmedt C, Wirapati P, Lallemand F, Tutt AM, Gillet C, 
#'   Ellis P, Ryder K, Reid JF, Daidone MG, Pierotti MA, Berns EMJJ, Jansen MPHM,
#'   Foekens JA, Delorenzi M, Bontempi G, Piccart MJ and Sotiriou C (2008) 
#'   "Predicting prognosis using molecular profiling in estrogen receptor-
#'   positive breast cancer treated with tamoxifen", BMC Genomics, 9(1):239
#'
#' @seealso
#' [genefu::gene76]
#' 
#' @examples
#' # load TAMR13 signature
#' data(sig.tamr13)
#' # load VDX dataset
#' data(vdxs)
#' # compute relapse score
#' tamr13.vdxs <- tamr13(data=data.vdxs, annot=annot.vdxs, do.mapping=FALSE)
#' summary(tamr13.vdxs$score)
#'
#' @md
#' @export
if(getRversion() >= "2.15.1")  utils::globalVariables("sig.tamr13")

tamr13 <-
function(data, annot, do.mapping=FALSE, mapping, verbose=FALSE) {

	rest <- resc <- NULL
	for(i in 1:length(sig.tamr13)) {
		rest <- cbind(rest, sig.score(x=sig.tamr13[[i]], data=data, annot=annot, do.mapping=do.mapping, mapping=mapping, signed=TRUE, verbose=verbose)$score)
		resc <- c(resc, attributes(sig.tamr13[[i]])$coefficient)
	}
	dimnames(rest) <- list(dimnames(data)[[1]], names(sig.tamr13))
	names(resc) <- names(sig.tamr13)
	res <- drop(rest %*% resc)
	riskt <- rep(NA, length(res))
	names(riskt) <- names(res)
	res <- list("score"=res, "risk"=riskt)
	
	return (res)
}
