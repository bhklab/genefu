#' @title Function to compute the PIK3CA gene signature (PIK3CA-GS)
#'
#' @description
#' This function computes signature scores from gene expression values
#'   following the algorithm used for the PIK3CA gene signature (PIK3CA-GS).
#'
#' @usage
#' pik3cags(data, annot, do.mapping = FALSE, mapping, verbose = FALSE)
#'
#' @param data Matrix of gene expressions with samples in rows and probes in
#'   columns, dimnames being properly defined.
#' @param annot Matrix of annotations with at least one column named
#'   "EntrezGene.ID", dimnames being properly defined.
#' @param do.mapping TRUE if the mapping through Entrez Gene ids must be
#'   performed (in case of ambiguities, the most variant probe is kept for
#'   each gene), FALSE otherwise.
#' @param mapping Matrix with columns "EntrezGene.ID" and "probe" used to force
#'   the mapping such that the probes are not selected based on their variance.
#' @param verbose TRUE to print informative messages, FALSE otherwise.
#'
#' @return
#' Vector of signature scores for PIK3CA-GS
#'
#' @references
#' Loi S, Haibe-Kains B, Majjaj S, Lallemand F, Durbecq V, Larsimont D,
#'   Gonzalez-Angulo AM, Pusztai L, Symmans FW, Bardelli A, Ellis P, Tutt AN,
#'   Gillett CE, Hennessy BT., Mills GB, Phillips WA, Piccart MJ, Speed TP,
#'   McArthur GA, Sotiriou C (2010) "PIK3CA mutations associated with gene
#'   signature of low mTORC1 signaling and better outcomes in estrogen
#'   receptor-positive breast cancer", Proceedings of the National Academy of
#'   Sciences, 107(22):10208-10213
#'
#' @seealso
#' [genefu::gene76]
#'
#' @examples
#' # load GGI signature
#' data(sig.pik3cags)
#' # load NKI dataset
#' data(nkis)
#' # compute relapse score
#' pik3cags.nkis <- pik3cags(data=data.nkis, annot=annot.nkis, do.mapping=TRUE)
#' head(pik3cags.nkis)
#'
#' @md
#' @export
#' @name pik3cags
pik3cags <- function(data, annot, do.mapping=FALSE, mapping, verbose=FALSE) {

	if (!exists('sig.pik3cags')) data(sig.pik3cags, envir=environment())

	pik3cags.gl <- sig.pik3cags[ ,c("probe", "EntrezGene.ID", "coefficient")]
	res <- sig.score(x=pik3cags.gl, data=data, annot=annot, do.mapping=do.mapping, mapping=mapping, signed=TRUE, verbose=verbose)$score

	return (res)
}
