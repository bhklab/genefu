#' @title Function to write a 'csv' file containing gene lists 
#'   (aka gene signatures)
#'
#' @description
#' This function allows for writing a 'csv' file containing gene signatures. 
#'   Each gene signature is composed of at least four columns: "gene.list" is 
#'   the name of the signature on the first line and empty fields below, 
#'   "probes" are the probe names, "EntrezGene.ID" are the EntrezGene IDs 
#'   and "coefficient" are the coefficients of each probe.
#'
#' @usage
#' write.m.file(obj, file, ...)
#'
#' @param obj	List of gene signatures.
#' @param file Filename of the 'csv' file.
#' @param ...	Additional parameters for read.csv function.
#'
#' @return
#' None.
#'
#' @examples
#' # load gene modules published by Demsedt et al 2009
#' data(mod1)
#' # write these gene modules in a 'csv' file
#' # Not run: write.m.file(obj=mod1, file="desmedt2009_genemodules.csv")
#' 
#' @md
#' @export
write.m.file <-
function(obj, file, ...) {
	lcn <- dimnames(obj[[1]])[[2]]
	c1 <- c2 <- NULL
	for (i in 1:length(obj)) {
		ct <- names(obj)[i]
		tt <- NULL
		for(j in 1:ncol(obj[[i]])) { tt <- cbind(tt, as.character(obj[[i]][ ,j]))}
		colnames(tt) <- colnames(obj[[i]])
		c1 <- c(c1, ct, rep("", nrow(tt)))
		c2 <- rbind(c2, tt, rep("", ncol(tt)))
	}
	dimnames(c2)[[1]] <- 1:nrow(c2)
	res <- cbind(c1, c2)[-length(c1), ,drop=FALSE]
	dimnames(res)[[2]] <- c("gene.list", lcn)
	dimnames(res)[[1]] <- 1:nrow(res)
	write.table(res, file=file, row.names=FALSE, sep=",", ...)
	invisible(res)
}