#' @title Function to compute the prediction strength of a clustering model
#'
#' @description
#' This function computes the prediction strength of a clustering model as published 
#'   in R. Tibshirani and G. Walther 2005.
#'
#' @usage
#' ps.cluster(cl.tr, cl.ts, na.rm = FALSE)
#'
#' @param cl.tr	Clusters membership as defined by the original clustering model, i.e. 
#'   the one that was not fitted on the dataset of interest.
#' @param cl.ts	Clusters membership as defined by the clustering model fitted on the 
#'   dataset of interest.
#' @param na.rm	TRUE if missing values should be removed, FALSE otherwise.
#'
#' @return
#' A list with items:
#' - ps: the overall prediction strength (minimum of the prediction strengths at cluster level).
#' - ps.cluster: Prediction strength for each cluster
#' - ps.individual: Prediction strength for each sample.
#'
#' @references
#' R. Tibshirani and G. Walther (2005) "Cluster Validation by Prediction Strength", 
#'   Journal of Computational and Graphical Statistics, 14(3):511-528.
#'
#' @examples
#' # load SSP signature published in Sorlie et al. 2003
#' data(ssp2003)
#' # load NKI data
#' data(nkis)
#' # SP2003 fitted on NKI
#' ssp2003.2nkis <- intrinsic.cluster(data=data.nkis, annot=annot.nkis,
#'   do.mapping=TRUE, std="robust",
#'   intrinsicg=ssp2003$centroids.map[ ,c("probe", "EntrezGene.ID")],
#'   number.cluster=5, mins=5, method.cor="spearman",
#'   method.centroids="mean", verbose=TRUE)
#' # SP2003 published in Sorlie et al 2003 and applied in VDX
#' ssp2003.nkis <- intrinsic.cluster.predict(sbt.model=ssp2003,
#'   data=data.nkis, annot=annot.nkis, do.mapping=TRUE, verbose=TRUE)
#' # prediction strength of sp2003 clustering model
#' ps.cluster(cl.tr=ssp2003.2nkis$subtype, cl.ts=ssp2003.nkis$subtype,
#'   na.rm = FALSE)
#'
#' @md
#' @export
ps.cluster <-
function(cl.tr, cl.ts, na.rm=FALSE) {
	## consider cl.ts as reference
	if(length(cl.tr) != length(cl.ts)) { stop("the two clustering must have the same length!") }
	if(is.null(names(cl.tr))) { names(cl.tr) <- names(cl.ts) <- paste("X", 1:length(cl.tr), sep=".") }
	cc.ix <- complete.cases(cl.tr, cl.ts)
	if(any(!cc.ix) & !na.rm) { stop("missing values are present!") }
	nn <- sum(cc.ix)
	cltr <- cl.tr[cc.ix]
	clts <- cl.ts[cc.ix]
	ucltr <- sort(unique(cltr))	
	uclts <- sort(unique(clts))
	if(length(ucltr) != length(uclts)) {
		message("the number of clusters should be the same")
		tt <- rep(NA, max(length(ucltr), length(uclts)))
		names(tt) <- ifelse(length(ucltr) > length(uclts), ucltr, uclts)
		tt2 <- rep(NA, length(cl.tr))
		names(tt2) <- names(cl.tr)
		return(list("ps"=0, "ps.cluster"=tt, "ps.individual"=tt2))
	}
	if(!all(ucltr == uclts)) { ## the number of clusters is the same but the labels differ
		message("the labels of the clusters differ")
		## keep labels from cl.ts
		tixtr <- !is.element(ucltr, uclts)
		tixts <- !is.element(uclts, ucltr)
		for(mm in 1:sum(tixtr)) {
			cltr[cltr == ucltr[which(tixtr)[mm]]] <- uclts[which(tixts)[mm]]
		}
		ucltr <- sort(unique(cltr))	
	}
	ll <- length(uclts)
	
	## co-membership matrix for the two clusterings
	Dtr <- matrix(NA, nrow=nn, ncol=nn, dimnames=list(names(cltr), names(cltr)))
	Dts <- matrix(NA, nrow=nn, ncol=nn, dimnames=list(names(clts), names(clts)))
	for(i in 1:(nn-1)) {
		for(j in (i+1):nn) {
			if(cltr[i] == cltr[j]) { Dtr[i,j] <- Dtr[j,i] <- 1 } else { Dtr[i,j] <- Dtr[j,i] <- 0 }
			if(clts[i] == clts[j]) { Dts[i,j] <- Dts[j,i] <- 1 } else { Dts[i,j] <- Dts[j,i] <- 0 }
		}
	}
	
	nltr <- table(cltr)[as.character(ucltr)]
	nlts <- table(clts)[as.character(uclts)]
	myps <- NULL
	for(l in 1:ll) {
		al <- which(clts == uclts[l])
		if(length(al) > 1) { rr <- (1 / (nlts[l] * (nlts[l] - 1))) * sum(Dtr[al,al], na.rm=TRUE) } else { rr <- 0}
		myps <- c(myps, rr)	
	}
	if(length(uclts) >= length(ucltr)) { names(myps) <- uclts } else { names(myps) <- ucltr }
	
	myps.ind <- NULL
	for(i in 1:nn) {
		aki <- Dts[i, ] == 1 & !is.na(Dts[i, ])
		rr <- (1 / sum(aki)) * sum(Dtr[i,aki] == 1, na.rm=TRUE)
		myps.ind <- c(myps.ind, rr)
	}
	names(myps.ind) <- names(clts)
	myps.ind2 <- rep(NA, length(cl.tr))
	names(myps.ind2) <- names(cl.tr)
	myps.ind2[names(myps.ind)] <- myps.ind

	return(list("ps"=min(myps), "ps.cluster"=myps, "ps.individual"=myps.ind2))
}
