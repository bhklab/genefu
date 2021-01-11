if (getRversion() >= "2.15.1")  utils::globalVariables("sig.oncotypedx")

#' @title Function to compute the OncotypeDX signature as published by
#'   Paik et al. in 2004.
#'
#' @description
#' This function computes signature scores and risk classifications from
#'   gene expression values following the algorithm used for the OncotypeDX
#'   signature as published by Paik et al. 2004.
#'
#' @usage
#' oncotypedx(data, annot, do.mapping = FALSE, mapping, verbose = FALSE)
#'
#' @param data Matrix of gene expressions with samples in rows and
#'   probes in columns, dimnames being properly defined.
#' @param annot Matrix of annotations with at least one column named
#'   "EntrezGene.ID", dimnames being properly defined.
#' @param do.mapping TRUE if the mapping through Entrez Gene ids must
#'   be performed (in case of ambiguities, the most variant probe is kept
#'   for each gene), FALSE otherwise. Note that for Affymetrix HGU
#'   datasets, the mapping is not necessary.
#' @param mapping Matrix with columns "EntrezGene.ID" and "probe" used
#'   to force the mapping such that the probes are not selected based on
#'   their variance.
#' @param verbose TRUE to print informative messages, FALSE otherwise.
#'
#' @details
#' Note that for Affymetrix HGU datasets, the mapping is not necessary.
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
#' S. Paik, S. Shak, G. Tang, C. Kim, J. Bakker, M. Cronin, F. L. Baehner,
#'   M. G. Walker, D. Watson, T. Park, W. Hiller, E. R. Fisher, D. L. Wickerham,
#'   J. Bryant, and N. Wolmark (2004) "A Multigene Assay to Predict Recurrence
#'   of Tamoxifen-Treated, Node-Negative Breast Cancer", New England Journal
#'   of Medicine, 351(27):2817-2826.
#'
#' @examples
#' # load GENE70 signature
#' data(sig.oncotypedx)
#' # load NKI dataset
#' data(nkis)
#' # compute relapse score
#' rs.nkis <- oncotypedx(data=data.nkis, annot=annot.nkis, do.mapping=TRUE)
#' table(rs.nkis$risk)
#'
#' @md
#' @export
oncotypedx <- function(data, annot, do.mapping=FALSE, mapping, do.scaling=TRUE,
                       verbose=FALSE) {

	## the reference genes are not taken into account due to their absence from most platforms
	sig2 <- sig.oncotypedx[sig.oncotypedx[ , "group"] != "reference",  , drop=FALSE]
	dimnames(sig2)[[1]] <- sig2[ , "probe.affy"]
	gt <- nrow(sig2)
	if(do.mapping) { ## not an affy HGU platform
		gid1 <- as.numeric(as.character(sig2[ ,"EntrezGene.ID"]))
		names(gid1) <- dimnames(sig2)[[1]]
		gid2 <- as.numeric(as.character(annot[ ,"EntrezGene.ID"]))
		names(gid2) <- dimnames(annot)[[1]]
		## remove missing and duplicated geneids from the gene list
		rm.ix <- is.na(gid1) | duplicated(gid1)
		gid1 <- gid1[!rm.ix]
		## mqpping
		rr <- geneid.map(geneid1=gid2, data1=data, geneid2=gid1, verbose=FALSE)
		gm <- length(rr$geneid2)
		mymapping <- c("mapped"=gm, "total"=gt)
		if(length(rr$geneid1) != gt) { ## some genes are missing
			res <- rep(NA, nrow(data))
			names(res) <- dimnames(data)[[1]]
			if(verbose) { message(sprintf("probe candidates: %i/%i", gm, gt)) }
			return(list("score"=res, "risk"=res, "mapping"=mymapping, "probe"=NA))
		}
		gid1 <- rr$geneid2
		gid2 <- rr$geneid1
		data <- rr$data1
		myprobe <- cbind("probe"=names(gid1), "EntrezGene.ID"=gid1, "new.probe"=names(gid2))
		## change the names of probes in the data
		dimnames(data)[[2]] <- names(gid2) <- names(gid1)
	} else {
		myprobe <- NA
		data <- data[ ,intersect(dimnames(sig2)[[1]], dimnames(data)[[2]])]
		gm <- ncol(data)
		mymapping <- c("mapped"=gm, "total"=gt)
		if(nrow(sig2) != ncol(data)) { ## some genes are missing
			res <- rep(NA, nrow(data))
			names(res) <- dimnames(data)[[1]]
			if(verbose) { message(sprintf("probe candidates: %i/%i", ncol(data), gt)) }
			return(list("score"=res, "risk"=res, "mapping"=mymapping, "probe"=myprobe))
		}
	}
	## rename gene names by the gene symbols
	dimnames(data)[[2]] <- dimnames(sig2)[[1]] <- sig2[ , "symbol"]

	if (do.scaling) {
		## scaling between 0 and 15
		data <- apply(data, 2, function(x) { xx <- (x - min(x, na.rm=TRUE)) / (max(x, na.rm=TRUE) - min(x, na.rm=TRUE)); return(xx * 15) })
	} else if (max(data, na.rm=TRUE) > 20 || min(data, na.rm=TRUE) < -5) {
		## check that each gene expression lies approximately in [0, 15]
		warning("The max and min values of your data suggest it is not already scaled...
			If it is please ignore this message, otherwise set `do.scaling=TRUE` to scale in the function call."
			)
	}

	## OcotypeDX recurrence score
	## GRB7 group score = 0.9 * GRB7 + 0.1 * HER2 if result < 8, then result = 8
	## ER group score = (0.8 * ER + 1.2 * PGR + BCL2 + SCUBE2) / 4
	## proliferation group score = ( survivin + KI67 + MYBL2 + CCNB1 + STK15) / 5 if result < 6.5, then result = 6.5
	## invasion group score = (CTSL2 + MMP11) / 2
	## RSU = + 0.47 * GRB7 group score - 0.34 * ER group score + 1.04 * proliferation group score + 0.10 * invasion group score + 0.05 * CD68 - 0.08 GSTM1 - 0.07 * BAG1

	cc.ix <- complete.cases(data)
	rs <- rs.unscaled <- rsrisk <- NULL
	for (i in 1:nrow(data)) {
		if(cc.ix[i]) {
			grb7.gs <- 0.9 * data[i, "GRB7"] + 0.1 * data[i, "ERBB2"]
			if (grb7.gs < 8) { grb7.gs <- 8 }

			er.gs <- (0.8 * data[i, "ESR1"] + 1.2 * data[i, "PGR"] + data[i, "BCL2"] + data[i, "SCUBE2"]) / 4

			proliferation.gs <- (data[i, "BIRC5"] + data[i, "MKI67"] + data[i, "MYBL2"] + data[i, "CCNB1"] + data[i, "AURKA"]) / 5
			if (proliferation.gs < 6.5) { proliferation.gs <- 6.5 }

			invasion.gs <- (data[i, "CTSL2"] + data[i, "MMP11"]) / 2

			rsu <- 0.47 * (grb7.gs) - 0.34 * (er.gs) + 1.04 * (proliferation.gs) + 0.1 * (invasion.gs) + 0.05 * data[i, "CD68"] - 0.08 * data[i, "GSTM1"] - 0.07 * data[i, "BAG1"]
			## rescale the score
			rsu2 <- rsu
			if(rsu >= 0 & rsu <= 100) { rsu <- 20 * (rsu - 6.7) }
			if(rsu < 0) { rsu <- 0 }
			if(rsu > 100) { rsu <- 100 }
			## use of the official curoffs
			if(rsu < 18) { rsr <- 0 }
			if(rsu >= 18 & rsu < 31) { rsr <- 0.5 }
			if(rsu >= 31) { rsr <- 1 }
		}
		else { rsu <- rsr <- rsu2 <- NA }
		rs.unscaled <- c(rs.unscaled, rsu2)
		rs <- c(rs, rsu)
		rsrisk <- c(rsrisk, rsr)
	}
	names(rs) <- names(rs.unscaled) <- names(rsrisk) <- dimnames(data)[[1]]
	return(list("score"=rs, "risk"=rsrisk, "mapping"=mymapping, "probe"=myprobe))
}
