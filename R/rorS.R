#' @title Function to compute the rorS signature as published by Parker
#'   et al 2009
#'
#' @description
#' This function computes signature scores and risk classifications from gene
#'   expression values following the algorithm used for the rorS signature as
#'   published by Parker et al 2009.
#'
#' @usage
#' rorS(data, annot, do.mapping = FALSE, mapping, verbose = FALSE)
#'
#' @param data	Matrix of gene expressions with samples in rows and
#'   probes in columns, dimnames being properly defined.
#' @param annot	Matrix of annotations with at least one column named
#'   "EntrezGene.ID", dimnames being properly defined.
#' @param do.mapping	TRUE if the mapping through Entrez Gene ids must be
#'   performed (in case of ambiguities, the most variant probe is kept for
#'   each gene), FALSE otherwise. Note that for Affymetrix HGU datasets, the
#'   mapping is not necessary.
#' @param mapping	Matrix with columns "EntrezGene.ID" and "probe" used to
#'   force the mapping such that the probes are not selected based on their
#'   variance.
#' @param verbose	TRUE to print informative messages, FALSE otherwis.
#'
#' @return
#' A list with items:
#' - score: Continuous signature scores
#' - risk: Binary risk classification, 1 being high risk and 0 being low risk.
#' - mapping: Mapping used if necessary.
#' - probe: If mapping is performed, this matrix contains the correspondence
#'   between the gene list (aka signature) and gene expression data.
#'
#' @references
#' Parker, Joel S. and Mullins, Michael and Cheang, Maggie C.U. and Leung,
#'   Samuel and Voduc, David and Vickery, Tammi and Davies, Sherri and Fauron,
#'   Christiane and He, Xiaping and Hu, Zhiyuan and Quackenbush, John F. and
#'   Stijleman, Inge J. and Palazzo, Juan and Marron, J.S. and Nobel,
#'   Andrew B. and Mardis, Elaine and Nielsen, Torsten O. and Ellis,
#'   Matthew J. and Perou, Charles M. and Bernard, Philip S. (2009) "Supervised
#'   Risk Predictor of Breast Cancer Based on Intrinsic Subtypes", Journal of
#'   Clinical Oncology, 27(8):1160-1167
#'
#' @examples
#' # load NKI dataset
#' data(vdxs)
#' data(pam50)
#'
#' # compute relapse score
#' rs.vdxs <- rorS(data=data.vdxs, annot=annot.vdxs, do.mapping=TRUE)
#'
#' @md
#' @export
#' @name rorS
rorS <- function(data, annot, do.mapping=FALSE, mapping, verbose=FALSE) {

  ## PAM50 classification
  data(pam50, envir=environment())
  sbts <- intrinsic.cluster.predict(sbt.model=pam50, data=data, annot=annot, do.mapping=do.mapping, verbose=FALSE)
  mymapping <- c("mapped"=nrow(sbts$centroids.map), "total"=nrow(pam50$centroids.map))
  ## ROR-S
  rs.unscaled <- rs <- rsrisk <- rep(NA, nrow(data))
  names(rs.unscaled) <- names(rs) <- names(rsrisk) <- rownames(data)
  rst <- 0.05 * sbts$cor[ , "Basal"] + 0.12 * sbts$cor[ , "Her2"] - 0.34 * sbts$cor[ , "LumA"] + 0.23 * sbts$cor[ , "LumB"]
  rs.unscaled[names(rst)] <- rst
  ## rescale between 0 and 100
  rs <- (rs.unscaled - quantile(rs.unscaled, probs=0.025, na.rm=TRUE)) / (quantile(rs.unscaled, probs=0.975, na.rm=TRUE) - quantile(rs.unscaled, probs=0.025, na.rm=TRUE)) * 100
  rs[!is.na(rs) & rs < 0] <- 0
  rs[!is.na(rs) & rs > 100] <- 100
  rsrisk[rs < 29] <- "Low"
  rsrisk[rs >= 29 & rs < 53] <- "Intermediate"
  rsrisk[rs >= 53] <- "High"
  rsrisk <- factor(rsrisk, levels=c("Low", "Intermediate", "High"))

	return(list("score"=rs, "risk"=rsrisk, "mapping"=mymapping, "probe"=sbts$centroids.map))
}
