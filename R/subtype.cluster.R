#' @title Function to fit the Subtype Clustering Model
#'
#' @description
#' This function fits the Subtype Clustering Model as published in Desmedt
#'   et al. 2008 and Wiarapati et al. 2008. This model is actually a mixture
#'   of three Gaussians with equal shape, volume and variance (see EEI model
#'   in Mclust). This model is adapted to breast cancer and uses ESR1, ERBB2
#'   and AURKA dimensions to identify the molecular subtypes, i.e. ER-/HER2-,
#'   HER2+ and ER+/HER2- (Low and High Prolif).
#'
#' @usage
#' subtype.cluster(module.ESR1, module.ERBB2, module.AURKA, data, annot,
#'   do.mapping = FALSE, mapping, do.scale = TRUE, rescale.q = 0.05,
#'   model.name = "EEI", do.BIC = FALSE, plot = FALSE, filen, verbose = FALSE)
#'
#' @param module.ESR1	Matrix containing the ESR1-related gene(s) in
#'   rows and at least three columns: "probe", "EntrezGene.ID" and
#'   "coefficient" standing for the name of the probe, the NCBI Entrez
#'   Gene id and the coefficient giving the direction and the strength of
#'   the association of each gene in the gene list.
#' @param module.ERBB2	Idem for ERBB2.
#' @param module.AURKA	Idem for AURKA.
#' @param data	Matrix of gene expressions with samples in rows and probes
#'   in columns, dimnames being properly defined.
#' @param annot	Matrix of annotations with at least one column named
#'   "EntrezGene.ID", dimnames being properly defined.
#' @param do.mapping	TRUE if the mapping through Entrez Gene ids must
#'   be performed (in case of ambiguities, the most variant probe is kept
#'   for each gene), FALSE otherwise.
#' @param mapping	**DEPRECATED** Matrix with columns "EntrezGene.ID" and
#'   "probe" used to force the mapping such that the probes are not selected
#'   based on their variance.
#' @param do.scale TRUE if the ESR1, ERBB2 and AURKA (module) scores must be
#'   rescaled (see rescale), FALSE otherwise.
#' @param rescale.q	Proportion of expected outliers for rescaling the gene expressions.
#' @param do.BIC	TRUE if the Bayesian Information Criterion must be computed
#'   for number of clusters ranging from 1 to 10, FALSE otherwise.
#' @param model.name Name of the model used to fit the mixture of Gaussians
#'   with the Mclust from the mclust package; default is "EEI" for fitting a
#'   mixture of Gaussians with diagonal variance, equal volume, equal shape
#'   and identical orientation.
#' @param plot TRUE if the patients and their corresponding subtypes must
#'   be plotted, FALSE otherwise.
#' @param filen Name of the csv file where the subtype clustering model must
#' be stored.
#' @param verbose TRUE to print informative messages, FALSE otherwise.
#'
#' @return
#' A list with items:
#' - model: Subtype Clustering Model (mixture of three Gaussians),
#'   like scmgene.robust, scmod1.robust and scmod2.robust when this function
#'   is used on expO dataset (International Genomics Consortium) with the gene
#'   modules published in the two references cited below.
#' - BIC: Bayesian Information Criterion for the Subtype Clustering Model
#'   with number of clusters ranging from 1 to 10.
#' - subtype: Subtypes identified by the Subtype Clustering Model. Subtypes
#'   can be either "ER-/HER2-", "HER2+" or "ER+/HER2-".
#' - subtype.proba: Probabilities to belong to each subtype estimated by
#'   the Subtype Clustering Model.
#' - subtype2: Subtypes identified by the Subtype Clustering Model using
#'   AURKA to discriminate low and high proliferative tumors. Subtypes can
#'   be either "ER-/HER2-", "HER2+", "ER+/HER2- High Prolif" or "ER+/HER2- Low Prolif".
#' - subtype.proba2: Probabilities to belong to each subtype (including
#'   discrimination between lowly and highly proliferative ER+/HER2- tumors,
#'   see subtype2) estimated by the Subtype Clustering Model.
#' - module.scores: Matrix containing ESR1, ERBB2 and AURKA module scores.
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
#'
#' @seealso
#' [genefu::subtype.cluster.predict()], [genefu::intrinsic.cluster()],
#' [genefu::intrinsic.cluster.predict()], [genefu::scmod1.robust],
#' [genefu::scmod2.robust]
#'
#' @examples
#' # example without gene mapping
#' # load expO data
#' data(expos)
#' # load gene modules
#' data(mod1)
#' # fit a Subtype Clustering Model
#' scmod1.expos <- subtype.cluster(module.ESR1=mod1$ESR1, module.ERBB2=mod1$ERBB2,
#'   module.AURKA=mod1$AURKA, data=data.expos, annot=annot.expos, do.mapping=FALSE,
#'   do.scale=TRUE, plot=TRUE, verbose=TRUE)
#' str(scmod1.expos, max.level=1)
#' table(scmod1.expos$subtype2)
#'
#' # example with gene mapping
#' # load NKI data
#' data(nkis)
#' # load gene modules
#' data(mod1)
#' # fit a Subtype Clustering Model
#' scmod1.nkis <- subtype.cluster(module.ESR1=mod1$ESR1, module.ERBB2=mod1$ERBB2,
#'   module.AURKA=mod1$AURKA, data=data.nkis, annot=annot.nkis, do.mapping=TRUE,
#'   do.scale=TRUE, plot=TRUE, verbose=TRUE)
#' str(scmod1.nkis, max.level=1)
#' table(scmod1.nkis$subtype2)
#'
#' @md
#' @import mclust
#' @export
subtype.cluster <- function(module.ESR1, module.ERBB2, module.AURKA, data,
	annot, do.mapping=FALSE, mapping, do.scale=TRUE, rescale.q=0.05,
	model.name="EEI", do.BIC=FALSE, plot=FALSE, filen, verbose=FALSE)
{
	if(missing(data) || missing(annot)) { stop("data, and annot parameters must be specified") }

	sbtn <- c("ER-/HER2-", "HER2+", "ER+/HER2-")
	sbtn2 <- c("ER-/HER2-", "HER2+", "ER+/HER2- High Prolif", "ER+/HER2- Low Prolif")
	sig.esr1 <- sig.score(x=module.ESR1, data=data, annot=annot, do.mapping=do.mapping, mapping=mapping, verbose=FALSE)
	sig.erbb2 <- sig.score(x=module.ERBB2, data=data, annot=annot, do.mapping=do.mapping, mapping=mapping, verbose=FALSE)
	sig.aurka <- sig.score(x=module.AURKA, data=data, annot=annot, do.mapping=do.mapping, mapping=mapping, verbose=FALSE)
	dd <- cbind("ESR1"=sig.esr1$score, "ERBB2"=sig.erbb2$score, "AURKA"=sig.aurka$score)
	rnn <- rownames(dd)
	m.mod <- list("ESR1"=cbind("probe"=as.character(sig.esr1$probe[ ,"new.probe"]), "EntrezGene.ID"=as.character(sig.esr1$probe[ ,"EntrezGene.ID"]), "coefficient"=module.ESR1[match(sig.esr1$probe[ ,"probe"], module.ESR1[ ,"probe"]), "coefficient"]), "ERBB2"=cbind("probe"=as.character(sig.erbb2$probe[ ,"new.probe"]), "EntrezGene.ID"=as.character(sig.erbb2$probe[ ,"EntrezGene.ID"]), "coefficient"=module.ERBB2[match(sig.erbb2$probe[ ,"probe"], module.ERBB2[ ,"probe"]), "coefficient"]), "AURKA"=cbind("probe"=as.character(sig.aurka$probe[ ,"new.probe"]), "EntrezGene.ID"=as.character(sig.aurka$probe[ ,"EntrezGene.ID"]), "coefficient"=module.AURKA[match(sig.aurka$probe[ ,"probe"], module.AURKA[ ,"probe"]), "coefficient"]))
	if(do.scale) {
		## the rescaling needs a large sample size!!!
		## necessary if we want to validate the classifier using a different dataset
		## the estimation of survival probabilities depends on the scale of the score
		dd <- apply(dd, 2, function(x) { return((rescale(x, q=rescale.q, na.rm=TRUE) - 0.5) * 2) })
		rownames(dd) <- rnn
	} else { rescale.q <- NA }
	rownames(dd) <- rownames(data)
	dd2 <- dd

	cc.ix <- complete.cases(dd[ , c("ESR1", "ERBB2"), drop=FALSE])
	dd <- dd[cc.ix, , drop=FALSE]
	if(all(!cc.ix)) { stop("None of ESR1 and ERBB2 genes are present!") }

	if(do.BIC) {
		## save the BIC values for the all the methods and a number of clusters from 1 to 10
		cluster.bic <- mclust::Mclust(data=dd[ , c("ESR1", "ERBB2"), drop=FALSE], modelNames=model.name, G=1:10)$BIC
	} else { cluster.bic <- NA }

	#identify the 3 subtypes
	rr3 <- mclust::Mclust(data=dd[ , c("ESR1", "ERBB2"), drop=FALSE], modelNames=model.name, G=3)
	#redefine classification to be coherent with subtypes
	uclass <- sort(unique(rr3$classification))
	uclass <- uclass[!is.na(uclass)]
	if(length(uclass) != 3) { stop("less than 3 subtypes are identified!") }
	mm <- NULL
	for(i in 1:length(uclass)) {
		mm <- c(mm, median(dd[rr3$classification == uclass[i],"ERBB2"], na.rm=TRUE) )
	}
	nclass <-  uclass[order(mm, decreasing=TRUE)[1]]
	mm <- NULL
	for(i in 1:length(uclass[-nclass])) {
		mm <- c(mm, median(dd[rr3$classification == uclass[-nclass][i],"ESR1"], na.rm=TRUE))
	}
	nclass <- c(uclass[-nclass][order(mm, decreasing=TRUE)[2]], nclass, uclass[-nclass][order(mm, decreasing=TRUE)[1]])
	#nclass contains the new order
	rr3$z <- rr3$z[ ,nclass, drop=FALSE]
	ncl <- rr3$classification
	for(i in 1:length(uclass)) {
		ncl[rr3$classification == nclass[i]] <- i
	}
	rr3$classification <- ncl
	rr3$parameters$pro <- rr3$parameters$pro[nclass]
	rr3$parameters$mean <- rr3$parameters$mean[ , nclass, drop=FALSE]
	rr3$parameters$variance$sigma <- rr3$parameters$variance$sigma[ , , nclass, drop=FALSE]

	if(plot) {
		if(do.scale) {
			myxlim <- myylim <- c(-2, 2)
		} else {
			myxlim <- range(dd[ , "ESR1"])
			myylim <- range(dd[ , "ERBB2"])
		}
		## plot the mixture of Gaussians of the model
		xx <- mclust:::grid1(50, range=myxlim)
		yy <- mclust:::grid1(50, range=myylim)
		xxyy <- mclust:::grid2(xx,yy)
		#density
		xyDens <- mclust::dens(modelName = rr3$modelName, data = xxyy, parameters = rr3$parameters)
		xyDens <- matrix(xyDens, nrow = length(xx), ncol = length(yy))
		par(pty = "s")
		zz <- xyDens
		#plot
		persp(x = xx, y = yy, z = zz, xlim=myxlim, ylim=myylim, theta=-25, phi=30, expand=0.5, xlab="ESR1", ylab="ERBB2", zlab="Density", col="darkgrey", ticktype="detailed")
	}

	## use the previously computed model to fit a new model in a supervised manner
	myclass <- mclust::unmap(rr3$classification)
	dimnames(myclass)[[1]] <- dimnames(dd)[[1]]
	mclust.tr <- mclust::mstep(modelName=model.name, data=dd[ , c("ESR1", "ERBB2"), drop=FALSE], z=myclass)
	dimnames(mclust.tr$z) <- dimnames(myclass)
	emclust.tr <- mclust::estep(modelName=model.name, data=dd[ , c("ESR1", "ERBB2"), drop=FALSE], parameters=mclust.tr$parameters)
	dimnames(emclust.tr$z) <- dimnames(myclass)
	class.tr <- mclust::map(emclust.tr$z, warn=FALSE)
	names(class.tr) <- dimnames(dd)[[1]]
	dimnames(mclust.tr$parameters$mean)[[2]] <- names(mclust.tr$parameters$pro) <- dimnames(mclust.tr$z)[[2]]

	## subtypes
	sbt <- rep(NA, nrow(data))
	names(sbt) <- dimnames(data)[[1]]
	sbt[names(class.tr)] <- sbtn[class.tr]
	sbt.proba <- matrix(NA, nrow(data), ncol=ncol(emclust.tr$z), dimnames=list(dimnames(data)[[1]], sbtn))
	sbt.proba[dimnames(emclust.tr$z)[[1]], ] <- emclust.tr$z
	## discriminate between luminal A and B using AURKA
	## since proliferation is a continuum we fit a Gaussian using AURKA expression of the ER+/HER2- tumors
	tt <- mclust::Mclust(dd2[complete.cases(sbt, dd2[ , "AURKA"]) & sbt == sbtn[3], "AURKA"], modelNames="E", G=1)
	gauss.prolif <- c("mean"=tt$parameters$mean, "sigma"=tt$parameters$variance$sigmasq)
	sbt2 <- sbt
	sbt2[sbt == sbtn[3]] <- NA
	## probability that tumor is highly proliferative
	pprolif <- pnorm(q=dd2[ , "AURKA"], mean=gauss.prolif["mean"], sd=gauss.prolif["sigma"], lower.tail=TRUE)
	## high proliferation
	sbt2[sbt == sbtn[3] & pprolif >= 0.5 & complete.cases(sbt, pprolif)] <- sbtn2[3]
	## low proliferation
	sbt2[sbt == sbtn[3] & pprolif < 0.5 & complete.cases(sbt, pprolif)] <- sbtn2[4]
	## subtype probabilities for luminal B and A
	sbt.proba2 <- matrix(NA, nrow(data), ncol=ncol(emclust.tr$z) + 1, dimnames=list(dimnames(data)[[1]], sbtn2))
	tt <- sbt.proba[ , sbtn[3]]
	tt2 <- pprolif
	tt <- cbind(tt * tt2, tt * (1 - tt2))
	colnames(tt) <- sbtn2[3:4]
	sbt.proba2[ , sbtn2[1:2]] <- sbt.proba[ , sbtn[1:2]]
	sbt.proba2[ , sbtn2[3:4]] <- tt[ , sbtn2[3:4]]

	if(plot) {
		## plot the clusters
		mclust::mclust2Dplot(data=dd[ , c("ESR1", "ERBB2"), drop=FALSE], what="classification", classification=class.tr, parameters=mclust.tr$parameters, colors=c("darkred", "darkgreen", "darkblue"), xlim=myxlim, ylim=myylim)
		legend(x="topleft", col=c("darkred", "darkgreen", "darkblue"), legend=sbtn, pch=mclust::mclust.options("classPlotSymbols")[1:length(uclass)], bty="n")
		## plot the clusters with luminals A and B
		mycol <- mypch <- sbt2
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
		## display the three circles representing the Gaussians
		for(kk in 1:3) { mclust::mvn2plot(mu=mclust.tr$parameters$mean[ ,kk], sigma=mclust.tr$parameters$variance$sigma[ , , kk]) }
	}

	if(!missing(filen)) {
		#save model parameters in a csv file for reuse
		write(x=sprintf("# Benjamin Haibe-Kains. All rights reserved."), file=paste(filen, "csv", sep="."))
		write(x=sprintf("# model.name: %s", model.name), append=TRUE, file=paste(filen, "csv", sep="."))
		mymean <- t(mclust.tr$parameters$mean)
		if(is.null(dimnames(mymean)[[1]])) { dimnames(mymean)[[1]] <- 1:nrow(mymean) }
		for(i in 1:nrow(mymean)) { write(x=sprintf("# mean.%s: %g %g", dimnames(mymean)[[1]][i], mymean[i,1], mymean[i,2]), append=TRUE, file=paste(filen, "csv", sep=".")) }
		mysigma <- diag(mysigma <- mclust.tr$parameters$variance$sigma[ , ,1])
		write(x=sprintf("# sigma: %g %g", mysigma[1], mysigma[2]), append=TRUE, file=paste(filen, "csv", sep="."))
		mypro <- mclust.tr$parameters$pro
		write(x=sprintf("# pro: %g %g %g", mypro[1], mypro[2], mypro[3]), append=TRUE, file=paste(filen, "csv", sep="."))
		myscale <- mclust.tr$parameters$variance$scale
		write(x=sprintf("# scale: %g", myscale), append=TRUE, file=paste(filen, "csv", sep="."))
		myshape <- mclust.tr$parameters$variance$shape
		write(x=sprintf("# shape: %g %g", myshape[1], myshape[2]), append=TRUE, file=paste(filen, "csv", sep="."))
		write(x=sprintf("# gaussian.AURKA.mean: %g", gauss.prolif[1]), append=TRUE, file=paste(filen, "csv", sep="."))
		write(x=sprintf("# gaussian.AURKA.sigma: %g", gauss.prolif[2]), append=TRUE, file=paste(filen, "csv", sep="."))
		write(x=sprintf("# rescale.q: %g", rescale.q), append=TRUE, file=paste(filen, "csv", sep="."))
		write(paste("\"", c("module", dimnames(m.mod[[1]])[[2]]), "\"", collapse=",", sep=""), sep="", append=TRUE, file=paste(filen, "csv", sep="."))
		write.m.file(m.mod, file=paste(filen, "csv", sep="."), col.names=FALSE, append=TRUE)
	}

	return(list("model"=c(mclust.tr["parameters"], list("gaussian.AURKA"=gauss.prolif), list("rescale.q"=rescale.q), list("mod"=m.mod)), "BIC"=cluster.bic, "subtype"=sbt, "subtype.proba"=sbt.proba, "subtype2"=sbt2, "subtype.proba2"=sbt.proba2,  "module.scores"=dd2))
}
