#' @title Function to read a 'csv' file containing gene lists (aka gene signatures)
#'
#' @description
#' This function allows for reading a 'csv' file containing gene signatures. 
#'   Each gene signature is composed of at least four columns: "gene.list" is the name 
#'   of the signature on the first line and empty fields below, "probes" are the probe 
#'   names, "EntrezGene.ID" are the EntrezGene IDs and "coefficient" are the coefficients 
#'   of each probe.
#'
#' @usage
#' read.m.file(file, ...)
#'
#' @param file	Filename of the 'csv' file.
#' @param ... Additional parameters for read.csv function.
#'
#' @return
#' List of gene signatures.
#'
#' @seealso
#' [genefu::mod1], [genefu::mod2], 'extdata/desmedt2008_genemodules.csv', 'extdata/haibekains2009_sig_genius.csv'
#'
#' @examples
#' # read the seven gene modules as published in Desmedt et al 2008
#' genemods <- read.m.file(system.file("extdata/desmedt2008_genemodules.csv",
#'   package = "genefu"))
#' str(genemods, max.level=1)
#' # read the three subtype signtaures from GENIUS
#' geniusm <- read.m.file(system.file("extdata/haibekains2009_sig_genius.csv",
#'   package = "genefu"))
#' str(geniusm, max.level=1)
#' 
#' @md
#' @export
read.m.file <-
function(file, ...) {
	obj.file <- read.csv(file, stringsAsFactors=FALSE, ...)
	if(sum(!is.na(obj.file[ ,1]) & obj.file[ ,1] != "") > 0) {
		ix.delim <- c(which(obj.file[ ,1] != "")[-1]-1, nrow(obj.file) + 1)
		ix.f <- ix.l <- 1
		groups <- NULL
		npp <- np <- NULL
		for (i in 1:length(ix.delim)) {
			ix.l <- ix.delim[i] - 1
			np <- c(np, as.character(obj.file[ix.f,1]))
			groups <- c(groups, rep(i, ix.l - ix.f + 1))
			npp <- rbind(npp, obj.file[ix.f:ix.l,2:ncol(obj.file)])
			ix.f <- ix.l + 2
		}
		ugroups <- unique(groups)
		obj <- NULL
		for (j in 1:length(ugroups)) {
			obj <- c(obj, list(npp[groups == ugroups[j], ]))
		}
	names(obj) <- np
	} else { obj <- list("module"=obj)}
	return(obj)
}