`subtype.cluster` <-
function(module.ESR1, module.ERBB2, module.AURKA, data, annot, do.mapping=FALSE, mapping, do.scale=TRUE, rescale.q=0.05, model.name="EEI", do.BIC=FALSE, plot=FALSE, filen, verbose=FALSE) {
	require(mclust)
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
	}
	
	dd2 <- dd
	
	cc.ix <- complete.cases(dd[ , c("ESR1", "ERBB2"), drop=FALSE])
	dd <- dd[cc.ix, , drop=FALSE]
	if(all(!cc.ix)) { stop("None of ESR1 and ERBB2 genes are present!") }

	if(do.BIC) {
		## save the BIC values for the all the methods and a number of clusters from 1 to 10
		cluster.bic <- Mclust(data=dd[ , c("ESR1", "ERBB2"), drop=FALSE], modelNames=model.name, G=1:10)$BIC
	} else { cluster.bic <- NA }

	#identify the 3 subtypes
	rr3 <- Mclust(data=dd[ , c("ESR1", "ERBB2"), drop=FALSE], modelNames=model.name, G=3)
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
		#plot the mixture of Gaussians of the model
		xx <- grid1(50, range=myxlim) 
		yy <- grid1(50, range=myylim)
		xxyy <- grid2(xx,yy)
		#density
		xyDens <- dens(modelName = rr3$modelName, data = xxyy, parameters = rr3$parameters)
		xyDens <- matrix(xyDens, nrow = length(xx), ncol = length(yy))
		par(pty = "s") 
		zz <- xyDens
		#plot
		persp(x = xx, y = yy, z = zz, xlim=myxlim, ylim=myylim, theta=-25, phi=30, expand=0.5, xlab="ESR1", ylab="ERBB2", zlab="Density", col="darkgrey", ticktype="detailed")
	}
	
	## use the previously computed model to fit a new model in a supervised manner
	myclass <- unmap(rr3$classification)
	dimnames(myclass)[[1]] <- dimnames(dd)[[1]]
	mclust.tr <- mstep(modelName=model.name, data=dd[ , c("ESR1", "ERBB2"), drop=FALSE], z=myclass)
	dimnames(mclust.tr$z) <- dimnames(myclass)
	emclust.tr <- estep(modelName=model.name, data=dd[ , c("ESR1", "ERBB2"), drop=FALSE], parameters=mclust.tr$parameters)
	dimnames(emclust.tr$z) <- dimnames(myclass)
	class.tr <- map(emclust.tr$z)
	names(class.tr) <- dimnames(dd)[[1]]
	dimnames(mclust.tr$parameters$mean)[[2]] <- names(mclust.tr$parameters$pro) <- dimnames(mclust.tr$z)[[2]]
	
	## subtypes
	sbt <- rep(NA, nrow(data))
	names(sbt) <- dimnames(data)[[1]]
	sbt[names(class.tr)] <- sbtn[class.tr]
	sbt.proba <- matrix(NA, nrow(data), ncol=ncol(emclust.tr$z), dimnames=list(dimnames(data)[[1]], sbtn))
	sbt.proba[dimnames(emclust.tr$z)[[1]], ] <- emclust.tr$z
	## discriminate between luminal A and B using AURKA
	sbt2 <- sbt
	sbt2[sbt == sbtn[3]] <- NA
	cc <- median(dd2[sbt == sbtn[3] & !is.na(sbt), "AURKA"], na.rm=TRUE)
	sbt2[sbt == sbtn[3] & dd2[ , "AURKA"] >= cc & complete.cases(sbt, dd2[ , "AURKA"])] <- sbtn2[3]
	sbt2[sbt == sbtn[3] & dd2[ , "AURKA"] < cc & complete.cases(sbt, dd2[ , "AURKA"])] <- sbtn2[4]
	
	if(plot) {
		## plot the clusters
		mclust2Dplot(data=dd[ , c("ESR1", "ERBB2"), drop=FALSE], what="classification", classification=class.tr, parameters=mclust.tr$parameters, colors=c("darkred", "darkgreen", "darkblue"), xlim=myxlim, ylim=myylim)
		gplots::smartlegend(x="left", y="top", col=c("darkred", "darkgreen", "darkblue"), legend=sbtn, pch=.Mclust$classPlotSymbols[1:length(uclass)], bty="n", bg="white")
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
		gplots::smartlegend(x="left", y="top", col=c("darkred", "darkgreen", "darkorange", "darkviolet"), legend=sbtn2, pch=c(17, 0, 10, 10), bty="n", bg="white")
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
		mycutoff.AURKA <- cc
		write(x=sprintf("# cutoff.AURKA: %g", cc), append=TRUE, file=paste(filen, "csv", sep="."))
		write(x=sprintf("# rescale.q: %g", rescale.q), append=TRUE, file=paste(filen, "csv", sep="."))
		write(paste("\"", c("module", dimnames(m.mod[[1]])[[2]]), "\"", collapse=",", sep=""), sep="", append=TRUE, file=paste(filen, "csv", sep="."))
		write.m.file(m.mod, file=paste(filen, "csv", sep="."), col.names=FALSE, append=TRUE)
	}
	
	return(list("model"=c(mclust.tr["parameters"], list("cutoff.AURKA"=cc), list("rescale.q"=rescale.q), list("mod"=m.mod)), "BIC"=cluster.bic, "subtype"=sbt, "subtype.proba"=sbt.proba, "subtype2"=sbt2, "module.scores"=dd2))
}