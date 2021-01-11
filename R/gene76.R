#' @title Function to compute the Relapse Score as published by Wang et al. 2005
#' #'
#' @description
#' This function computes signature scores and risk classifications from gene
#'   expression values following the algorithm used for the Relapse Score (GENE76) as
#'   published by Wang et al. 2005.
#'
#' @usage
#' gene76(data, er)
#'
#' @param data Matrix of gene expressions with samples in rows and probes in columns,
#'   dimnames being properly defined.
#' @param er Vector containing the estrogen receptor (ER) status of breast cancer patients in
#'   the dataset.
#'
#'
#' @return
#' A list with items:
#' - score Continuous signature scores
#' - risk Binary risk classification, 1 being high risk and 0 being low risk.
#'
#' @references
#' Y. Wang and J. G. Klijn and Y. Zhang and A. M. Sieuwerts and M. P. Look and F.
#'   Yang and D. Talantov and M. Timmermans and M. E. Meijer-van Gelder and J. Yu and T.
#'   Jatkoe and E. M. Berns and D. Atkins and J. A. Foekens (2005) "Gene-Expression
#'   Profiles to Predict Distant Metastasis of Lymph-Node-Negative Primary Breast Cancer",
#'   Lancet, 365(9460):671–679.
#'
#' @seealso
#' [genefu::ggi]
#'
#' @examples
#' # load GENE76 signature
#' data(sig.gene76)
#' # load VDX dataset
#' data(vdxs)
#' # compute relapse score
#' rs.vdxs <- gene76(data=data.vdxs, er=demo.vdxs[ ,"er"])
#' table(rs.vdxs$risk)
#'
#' @md
#' @export
gene76 <- function(data, er) {

	A <- 313.5
	B <- 280

	score <- NULL
	for(i in 1:nrow(data)) {
		if(is.na(er[i])) { score <- c(score, NA) }
		else {
			if(er[i] == 1) {
				score <- c(score, A + sum(data[i,dimnames(sig.gene76)[[1]][sig.gene76[ ,"er"] == 1]] * sig.gene76[sig.gene76[ ,"er"] == 1,"std.cox.coefficient"]))
			}
			else {
				score <- c(score, B + sum(data[i,dimnames(sig.gene76)[[1]][sig.gene76[ ,"er"] == 0]] * sig.gene76[sig.gene76[ ,"er"] == 0,"std.cox.coefficient"]))
			}
		}
	}
	names(score) <- dimnames(data)[[1]]
	risk <- ifelse(score >= 0, 1, 0)

	return(list("score"=score, "risk"=risk))
}
