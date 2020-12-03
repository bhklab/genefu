#' @title Function to identify breast cancer molecular subtypes using 
#'   the Subtype Clustering Model
#'
#' @description
#' This function identifies the breast cancer molecular subtypes using a 
#'   Subtype Clustering Model fitted by subtype.cluster.
#'
#' @usage
#' subtype.cluster.predict(sbt.model, data, annot, do.mapping = FALSE,
#'   mapping, do.prediction.strength = FALSE,
#'   do.BIC = FALSE, plot = FALSE, verbose = FALSE)
#'
#' @param sbt.model	Subtype Clustering Model as returned by subtype.cluster.
#' @param data Matrix of gene expressions with samples in rows and probes in 
#'   columns, dimnames being properly defined.
#' @param annot	Matrix of annotations with at least one column named 
#'   "EntrezGene.ID", dimnames being properly defined.
#' @param do.mapping TRUE if the mapping through Entrez Gene ids must be 
#'   performed (in case of ambiguities, the most variant probe is kept 
#'   for each gene), FALSE otherwise.
#' @param mapping	**DEPRECATED** Matrix with columns "EntrezGene.ID" and 
#'   "probe" used to force the mapping such that the probes are not selected 
#'   based on their variance.
#' @param do.prediction.strength TRUE if the prediction strength must be 
#'   computed (Tibshirani and Walther 2005), FALSE otherwise.
#' @param do.BIC TRUE if the Bayesian Information Criterion must be computed 
#'   for number of clusters ranging from 1 to 10, FALSE otherwise.
#' @param plot TRUE if the patients and their corresponding subtypes must 
#'   be plotted, FALSE otherwise.
#' @param verbose	TRUE to print informative messages, FALSE otherwise.
#'
#' @return
#' A list with items:
#' - subtype: Subtypes identified by the Subtype Clustering Model. 
#'   Subtypes can be either "ER-/HER2-", "HER2+" or "ER+/HER2-".
#' - subtype.proba: Probabilities to belong to each subtype estimated 
#'   by the Subtype Clustering Model.
#' - prediction.strength: Prediction strength for subtypes.
#' - BIC: Bayesian Information Criterion for the Subtype Clustering Model 
#'   with number of clusters ranging from 1 to 10.
#' - subtype2: Subtypes identified by the Subtype Clustering Model using 
#'   AURKA to discriminate low and high proliferative tumors. Subtypes can be 
#'   either "ER-/HER2-", "HER2+", "ER+/HER2- High Prolif" or 
#'   "ER+/HER2- Low Prolif".
#' - subtype.proba2: Probabilities to belong to each subtype (including 
#'   discrimination between lowly and highly proliferative ER+/HER2- tumors, 
#'   see subtype2) estimated by the Subtype Clustering Model.
#' - prediction.strength2: Prediction strength for subtypes2.
#' - module.scores: Matrix containing ESR1, ERBB2 and AURKA module scores.
#' - mapping: Mapping if necessary (list of matrices with 3 columns: probe, 
#'   EntrezGene.ID and new.probe).
#' 
#' @references
#' Desmedt C, Haibe-Kains B, Wirapati P, Buyse M, Larsimont D, Bontempi G, 
#'   Delorenzi M, Piccart M, and Sotiriou C (2008) "Biological processes 
#'   associated with breast cancer clinical outcome depend on the molecular 
#'   subtypes", Clinical Cancer Research, 14(16):5158-5165.
#' Wirapati P, Sotiriou C, Kunkel S, Farmer P, Pradervand S, Haibe-Kains B, 
#'   Desmedt C, Ignatiadis M, Sengstag T, Schutz F, Goldstein DR, Piccart MJ 
#'   and Delorenzi M (2008) "Meta-analysis of Gene-Expression Profiles in 
#'   Breast Cancer: Toward a Unified Understanding of Breast Cancer Sub-typing 
#'   and Prognosis Signatures", Breast Cancer Research, 10(4):R65.
#' Tibshirani R and Walther G (2005) "Cluster Validation by Prediction 
#'   Strength", Journal of Computational and Graphical Statistics, 
#'   14(3):511-528
#'
#' @seealso
#' [genefu::subtype.cluster], [genefu::scmod1.robust], [genefu::scmod2.robust]
#' 
#' @examples
#' # without mapping (affy hgu133a or plus2 only)
#' # load VDX data
#' data(vdxs)
#' # Subtype Clustering Model fitted on EXPO and applied on VDX
#' sbt.vdxs <- subtype.cluster.predict(sbt.model=scmgene.robust, data=data.vdxs, 
#'   annot=annot.vdxs, do.mapping=FALSE, do.prediction.strength=FALSE, 
#'   do.BIC=FALSE, plot=TRUE, verbose=TRUE)
#' table(sbt.vdxs$subtype)
#' table(sbt.vdxs$subtype2)
#' 
#' # with mapping
#' # load NKI data
#' data(nkis)
#' # Subtype Clustering Model fitted on EXPO and applied on NKI
#' sbt.nkis <- subtype.cluster.predict(sbt.model=scmgene.robust, data=data.nkis, 
#'   annot=annot.nkis, do.mapping=TRUE, do.prediction.strength=FALSE, 
#'   do.BIC=FALSE, plot=TRUE, verbose=TRUE)
#' table(sbt.nkis$subtype)
#' table(sbt.nkis$subtype2)
#'
#' @md
#' @export
subtype.cluster.predict <-
function(sbt.model, data, annot, do.mapping=FALSE, mapping, do.prediction.strength=FALSE, do.BIC=FALSE, plot=FALSE, verbose=FALSE) {
	#require(mclust)
	if(missing(data) || missing(annot)) { stop("data, and annot parameters must be specified") }
	
	sbtn <- c("ER-/HER2-", "HER2+", "ER+/HER2-")
	sbtn2 <- c("ER-/HER2-", "HER2+", "ER+/HER2- High Prolif", "ER+/HER2- Low Prolif")
	
	if(is.list(sbt.model)) {
		## retrieve model
		subtype.c <- sbt.model[!is.element(names(sbt.model), c("cutoff.AURKA", "mod"))]
		model.name <- subtype.c$parameters$variance$modelName
		cc <- sbt.model$gaussian.AURKA
		mq <- sbt.model$rescale.q
		m.mod <- sbt.model$mod
	} else {
		## read model file
		rr <- readLines(con=sbt.model, n=11)[-1]
		nn <- unlist(lapply(X=rr, FUN=function(x) { x <- unlist(strsplit(x=unlist(strsplit(x=x, split=":"))[1], split=" ")); x <- x[length(x)]; return(x); }))
		m.param <- c(list(rr[[1]]), lapply(X=rr[-1], FUN=function(x) { x <- as.numeric(unlist(strsplit(x=unlist(strsplit(x=x, split=":"))[2], split=" "))); x <- x[!is.na(x)]; return(x)}))
		names(m.param) <- nn
		cn <- unlist(lapply(strsplit(nn[grep(pattern="mean", x=nn)], split="[.]"), FUN=function(x) { return(x[[2]]) }))
		m.mod <- read.m.file(sbt.model, comment.char="#")
		#construct a fake mclust object with the parameters of the model
		subtype.c <- NULL
		tt <- m.param$pro
		names(tt) <- cn
		subtype.c$parameters$pro <- tt
		tt <- sapply(X=m.param[grep(pattern="mean", x=nn)], FUN=function(x) { return(x) })
		dimnames(tt) <- list(names(m.mod)[1:2], cn)
		subtype.c$parameters$mean <- tt
		subtype.c$parameters$variance$modelName <- model.name <- m.param$modelname
		subtype.c$parameters$variance$d <- 2
		subtype.c$parameters$variance$G <- 3
		tt <- matrix(0, ncol=2, nrow=2, dimnames=list(names(m.mod)[1:2], names(m.mod)[1:2]))
		diag(tt) <- m.param$sigma
		subtype.c$parameters$variance$sigma <- array(tt, dim=c(2,2,3), dimnames=list(names(m.mod)[1:2], names(m.mod)[1:2], cn))
		subtype.c$parameters$variance$Sigma <- tt
		subtype.c$parameters$variance$scale <- m.param$scale
		subtype.c$parameters$variance$shape <- m.param$shape
		cc <- c("mean"=m.param$gaussian.AURKA.mean, "sigma"=m.param$gaussian.AURKA.sigma)
		mq <- m.param$rescale.q
	}
	do.scale <- ifelse(is.na(mq), FALSE, TRUE)
	sbt <- rep(NA, nrow(data))
	names(sbt) <- dimnames(data)[[1]]
	sbt.proba <- matrix(NA, nrow(data), ncol=length(sbtn), dimnames=list(dimnames(data)[[1]], sbtn))
	
	sigs.esr1 <- sig.score(x=m.mod$ESR1, data=data, annot=annot, do.mapping=do.mapping, mapping=mapping, verbose=FALSE)
	sigs.erbb2 <- sig.score(x=m.mod$ERBB2, data=data, annot=annot, do.mapping=do.mapping, mapping=mapping, verbose=FALSE)
	sigs.aurka <- sig.score(x=m.mod$AURKA, data=data, annot=annot, do.mapping=do.mapping, mapping=mapping, verbose=FALSE)
	## signature scores
	dd <- cbind("ESR1"=sigs.esr1$score, "ERBB2"=sigs.erbb2$score, "AURKA"=sigs.aurka$score)
	## mapping
	mymap <- list("ESR1"=sigs.esr1$probe, "ERBB2"=sigs.erbb2$probe, "AURLA"=sigs.aurka$probe)
	cln <- dimnames(subtype.c$parameters$mean)[[2]] <- as.character(1:ncol(subtype.c$parameters$mean))
	
	if(do.scale) {
		## the rescaling needs a large sample size!!!
		## necessary if we want to validate the classifier using a different dataset
		## the estimation of survival probabilities depends on the scale of the score
		dd <- apply(dd, 2, function(x) { return((rescale(x, q=mq, na.rm=TRUE) - 0.5) * 2) })
	}
	rownames(dd) <- rownames(data)
	dd2 <- dd
	
	cc.ix <- complete.cases(dd[ , c("ESR1", "ERBB2"), drop=FALSE])
	if(all(!cc.ix)) {
		ps.res <- ps.res2 <- BIC.res <- NULL
		if(do.prediction.strength) {
			tt <- rep(NA, length(sbtn))
			names(tt) <- sbtn
			tt2 <- rep(NA, length(sbtn2))
			names(tt2) <- sbtn2
			ps.res <- list("ps"=NA, "ps.cluster"=tt, "ps.individual"=sbt)
			ps.res2 <- list("ps"=NA, "ps.cluster"=tt2, "ps.individual"=sbt)
		}
		if(do.BIC) {
			BIC.res <- rep(NA, 10)
			names(BIC.res) <- 1:10
		}
		return(list("subtype"=sbt, "subtype.proba"=sbt.proba, "prediction.strength"=ps.res, "BIC"=BIC.res, "subtype2"=sbt, "prediction.strength2"=ps.res2))
	}
	dd <- dd[cc.ix, , drop=FALSE]
	
	emclust.ts <- mclust::estep(modelName=model.name, data=dd[ , c("ESR1", "ERBB2"), drop=FALSE], parameters=subtype.c$parameters)
	dimnames(emclust.ts$z) <- list(dimnames(dd)[[1]], cln)
	class.ts <- mclust::map(emclust.ts$z, warn=FALSE)
	names(class.ts) <- dimnames(dd)[[1]]
	uclass <- sort(unique(class.ts))
	uclass <- uclass[!is.na(uclass)]
	
	ps.res <- ps.res2 <- NULL
	if(do.prediction.strength) {
		if(nrow(dd) < 10) {
			warning("at least 10 observations are required to compute the prediction strength!")
			tt <- rep(NA, length(sbtn))
			names(tt) <- sbtn
			tt2 <- rep(NA, nrow(dd2))
			names(tt2) <- dimnames(dd2)[[1]]
			ps.res <- list("ps"=0, "ps.cluster"=tt, "ps.individual"=tt2)
			tt <- rep(NA, length(sbtn2))
			names(tt) <- sbtn2
			ps.res2 <- list("ps"=0, "ps.cluster"=tt, "ps.individual"=tt2)
		} else {
			## computation of the prediction strength of the clustering
			rr3 <- mclust::Mclust(data=dd[ , c("ESR1", "ERBB2"), drop=FALSE], modelNames=model.name, G=3)
			## redefine classification to be coherent with subtypes
			uclass <- sort(unique(rr3$classification))
			uclass <- uclass[!is.na(uclass)]
			if(length(uclass) != 3) {
				warning("less than 3 subtypes are identified!")
				tt <- rep(NA, length(sbtn))
				names(tt) <- sbtn
				tt2 <- rep(NA, nrow(dd2))
				names(tt2) <- dimnames(dd2)[[1]]
				ps.res <- list("ps"=0, "ps.cluster"=tt, "ps.individual"=tt2)
				tt <- rep(NA, length(sbtn2))
				names(tt) <- sbtn2
				ps.res2 <- list("ps"=0, "ps.cluster"=tt, "ps.individual"=tt2)
			} else {
				mm <- NULL
				for(i in 1:length(uclass)) {
					mm <- c(mm, median(dd[rr3$classification == uclass[i],"ERBB2"], na.rm=TRUE) )
				}
				nclass <-  uclass[order(mm, decreasing=TRUE)[1]]
				mm <- NULL
				for(i in 1:length(uclass[-nclass])) {
					mm <- c(mm, median(dd[rr3$classification == uclass[-nclass][i],"ESR1"], na.rm=TRUE) )
				}
				nclass <- c(uclass[-nclass][order(mm, decreasing=TRUE)[2]], nclass, uclass[-nclass][order(mm, decreasing=TRUE)[1]])
				## nclass contains the new order
				ncl <- rr3$classification
				for(i in 1:length(uclass)) {
					ncl[rr3$classification == nclass[i]] <- i
				}
				## use the previously computed model to fit a new model in a supervised manner
				myclass <- mclust::unmap(ncl)
				dimnames(myclass) <-  list(dimnames(dd)[[1]], sbtn)
				mclust.tr <- mclust::mstep(modelName=model.name, data=dd[ , c("ESR1", "ERBB2"), drop=FALSE], z=myclass)
				dimnames(mclust.tr$z) <- dimnames(myclass)
				emclust.tr <- mclust::estep(modelName=model.name, data=dd[ , c("ESR1", "ERBB2"), drop=FALSE], parameters=mclust.tr$parameters)
				dimnames(emclust.tr$z) <- dimnames(myclass)
				class.tr <- mclust::map(emclust.tr$z, warn=FALSE)
				names(class.tr) <- dimnames(dd)[[1]]
				## prediction strength
				ps.res <- ps.cluster(cl.tr=class.ts, cl.ts=class.tr, na.rm=TRUE)
				names(ps.res$ps.cluster) <- sbtn
				## check for missing values in ps.individual
				tt2 <- rep(NA, nrow(dd2))
				names(tt2) <- dimnames(dd2)[[1]]
				tt2[names(ps.res$ps.individual)] <- ps.res$ps.individual
				ps.res$ps.individual <- tt2
				## prediction strength with the separation in high and low proliferative tumors
				## since proliferation is a continuum we fit a Gaussian using AURKA expression of the ER+/HER2- tumors
				## refitted model
				tt <- mclust::Mclust(dd[complete.cases(class.tr, dd[ , "AURKA"]) & class.tr == 3, "AURKA"], modelNames="E", G=1)
				gauss.prolif <- c("mean"=tt$parameters$mean, "sigma"=tt$parameters$variance$sigmasq)
				class.tr2 <- class.tr
				class.tr2[class.tr == 3] <- NA
				## probability that tumor is highly proliferative
				pprolif <- pnorm(q=dd[ , "AURKA"], mean=gauss.prolif["mean"], sd=gauss.prolif["sigma"], lower.tail=TRUE)
				## high proliferation
				class.tr2[class.tr == 3 & pprolif >= 0.5 & complete.cases(class.tr, pprolif)] <- 3
				## low proliferation
				class.tr2[class.tr == 3 & pprolif < 0.5 & complete.cases(class.tr, pprolif)] <- 4
				## existing model
				tt <- mclust::Mclust(dd[complete.cases(class.ts, dd[ , "AURKA"]) & class.ts == 3, "AURKA"], modelNames="E", G=1)
				gauss.prolif <- c("mean"=tt$parameters$mean, "sigma"=tt$parameters$variance$sigmasq)
				class.ts2 <- class.ts
				class.ts2[class.ts == 3] <- NA
				## probability that tumor is highly proliferative
				pprolif <- pnorm(q=dd[ , "AURKA"], mean=gauss.prolif["mean"], sd=gauss.prolif["sigma"], lower.tail=TRUE)
				## high proliferation
				class.ts2[class.ts == 3 & pprolif >= 0.5 & complete.cases(class.ts, pprolif)] <- 3
				## low proliferation
				class.ts2[class.ts == 3 & pprolif < 0.5 & complete.cases(class.ts, pprolif)] <- 4
				## compute the prediction strength
				ps.res2 <- ps.cluster(cl.tr=class.ts2, cl.ts=class.tr2, na.rm=TRUE)
				names(ps.res2$ps.cluster) <- sbtn2
				## check for missing values in ps.individual
				tt2 <- rep(NA, nrow(dd2))
				names(tt2) <- dimnames(dd2)[[1]]
				tt2[names(ps.res2$ps.individual)] <- ps.res2$ps.individual
				ps.res2$ps.individual <- tt2
			}
		}
	}
	
	BIC.res <- NULL
	if(do.BIC) { 
		if(nrow(dd) >= 10) { BIC.res <- mclust::mclustBIC(data=dd[ , c("ESR1", "ERBB2"), drop=FALSE], modelNames=c(model.name), G=1:10)[ ,model.name] } else { warning("at least 10 observations are required to compute the BIC!") }
	}

	## subtypes
	sbt[names(class.ts)] <- sbtn[class.ts]
	sbt.proba[dimnames(emclust.ts$z)[[1]], ] <- emclust.ts$z
	## discriminate between luminal A and B using AURKA
	gauss.prolif <- cc
	sbt2 <- sbt
	sbt2[sbt == sbtn[3]] <- NA
	## probability that tumor is highly proliferative
	pprolif <- pnorm(q=dd2[ , "AURKA"], mean=gauss.prolif["mean"], sd=gauss.prolif["sigma"], lower.tail=TRUE)
	## high proliferation
	sbt2[sbt == sbtn[3] & pprolif >= 0.5 & complete.cases(sbt, pprolif)] <- sbtn2[3]
	## low proliferation
	sbt2[sbt == sbtn[3] & pprolif < 0.5 & complete.cases(sbt, pprolif)] <- sbtn2[4]
	## subtype probabilities for luminal B and A
	sbt.proba2 <- matrix(NA, nrow(data), ncol=length(sbtn2), dimnames=list(dimnames(data)[[1]], sbtn2))
	tt <- sbt.proba[ , sbtn[3]]
	tt2 <- pprolif
	tt <- cbind(tt * tt2, tt * (1 - tt2))
	colnames(tt) <- sbtn2[3:4]
	sbt.proba2[ , sbtn2[1:2]] <- sbt.proba[ , sbtn[1:2]]
	sbt.proba2[ , sbtn2[3:4]] <- tt[ , sbtn2[3:4]]
	
	if(plot) {
		if(do.scale) {
			myxlim <- myylim <- c(-2, 2)
		} else {
			myxlim <- range(dd[ , "ESR1"])
			myylim <- range(dd[ , "ERBB2"])
		}
		## plot the clusters with proliferation
		mycol <- mypch <- rep(NA, length(sbt2))
		mycol[sbt2 == sbtn2[1]] <- "darkred"
		mycol[sbt2 == sbtn2[2]] <- "darkgreen"
		mycol[sbt2 == sbtn2[3]] <- "darkorange"
		mycol[sbt2 == sbtn2[4]] <- "darkviolet"
		mypch[sbt2 == sbtn2[1]] <- 17
		mypch[sbt2 == sbtn2[2]] <- 0
		mypch[sbt2 == sbtn2[3] | sbt2 == sbtn2[4]] <- 10
		mypch <- as.numeric(mypch)
		names(mycol) <- names(mypch) <- names(sbt2)
		plot(x=dd[ , "ESR1"], y=dd[ , "ERBB2"], xlim=myxlim, ylim=myylim, xlab="ESR1", ylab="ERBB2", col=mycol[dimnames(dd)[[1]]], pch=mypch[dimnames(dd)[[1]]])
		legend(x="topleft", col=c("darkred", "darkgreen", "darkorange", "darkviolet"), legend=sbtn2, pch=c(17, 0, 10, 10), bty="n")
	}

	return(list("subtype"=sbt, "subtype.proba"=sbt.proba, "prediction.strength"=ps.res, "BIC"=BIC.res, "subtype2"=sbt2, "subtype.proba2"=sbt.proba2, "prediction.strength2"=ps.res2, "module.scores"=dd2, "mapping"=mymap))
}
