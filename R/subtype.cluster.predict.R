`subtype.cluster.predict` <-
function(sbt.model, data, annot, do.mapping=FALSE, mapping, do.scale=TRUE, do.prediction.strength=FALSE, do.BIC=FALSE, plot=FALSE, verbose=FALSE) {
	require(mclust)
	if(missing(data) || missing(annot)) { stop("data, and annot parameters must be specified") }
	
	sbtn <- c("ER-/HER2-", "HER2+", "ER+/HER2-")
	sbtn2 <- c("ER-/HER2-", "HER2+", "ER+/HER2- High Prolif", "ER+/HER2- Low Prolif")
	
	if(is.list(sbt.model)) {
		## retrieve model
		subtype.c <- sbt.model[!is.element(names(sbt.model), c("cutoff.AURKA", "mod"))]
		model.name <- subtype.c$parameters$variance$modelName
		cc <- sbt.model$cutoff.AURKA
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
		cc <- m.param$cutoff.AURKA
		mq <- m.param$rescale.q
	}
	
	sbt <- rep(NA, nrow(data))
	names(sbt) <- dimnames(data)[[1]]
	sbt.proba <- matrix(NA, nrow(data), ncol=3, dimnames=list(dimnames(data)[[1]], sbtn))
	
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
	
	emclust.ts <- estep(modelName=model.name, data=dd[ , c("ESR1", "ERBB2"), drop=FALSE], parameters=subtype.c$parameters)
	dimnames(emclust.ts$z) <- list(dimnames(dd)[[1]], cln)
	class.ts <- map(emclust.ts$z)
	names(class.ts) <- dimnames(dd)[[1]]
	uclass <- sort(unique(class.ts))
	uclass <- uclass[!is.na(uclass)]
	
	ps.res <- ps.res2 <- NULL
	if(do.prediction.strength) {
		#computation of the prediction strength of the clustering
		rr3 <- Mclust(data=dd[ , c("ESR1", "ERBB2"), drop=FALSE], modelNames=model.name, G=3)
		#redefine classification to be coherent with subtypes
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
			#nclass contains the new order
			ncl <- rr3$classification
			for(i in 1:length(uclass)) {
				ncl[rr3$classification == nclass[i]] <- i
			}
			#use the previously computed model to fit a new model in a supervised manner
			myclass <- unmap(ncl)
			dimnames(myclass) <-  list(dimnames(dd)[[1]], sbtn)
			mclust.tr <- mstep(modelName=model.name, data=dd[ , c("ESR1", "ERBB2"), drop=FALSE], z=myclass)
			dimnames(mclust.tr$z) <- dimnames(myclass)
			emclust.tr <- estep(modelName=model.name, data=dd[ , c("ESR1", "ERBB2"), drop=FALSE], parameters=mclust.tr$parameters)
			dimnames(emclust.tr$z) <- dimnames(myclass)
			class.tr <- map(emclust.tr$z)
			names(class.tr) <- dimnames(dd)[[1]]
			#prediction strength
			ps.res <- ps.cluster(cl.tr=class.ts, cl.ts=class.tr, na.rm=TRUE)
			names(ps.res$ps.cluster) <- sbtn
			## check for missing values in ps.individual
			tt2 <- rep(NA, nrow(dd2))
			names(tt2) <- dimnames(dd2)[[1]]
			tt2[names(ps.res$ps.individual)] <- ps.res$ps.individual
			ps.res$ps.individual <- tt2
			#prediction strength with the separation in high and low proliferative tumors
			class.ts2 <- class.ts
			class.ts2[class.ts == 3 & dd2[dimnames(dd)[[1]], "AURKA"] < cc & complete.cases(class.ts, dd2[dimnames(dd)[[1]], "AURKA"])] <- 4
			class.tr2 <- class.tr
			class.tr2[class.tr == 3 & dd2[dimnames(dd)[[1]], "AURKA"] < cc & complete.cases(class.tr, dd2[dimnames(dd)[[1]], "AURKA"])] <- 4
			ps.res2 <- ps.cluster(cl.tr=class.ts2, cl.ts=class.tr2, na.rm=TRUE)
			names(ps.res2$ps.cluster) <- sbtn2
			## check for missing values in ps.individual
			tt2 <- rep(NA, nrow(dd2))
			names(tt2) <- dimnames(dd2)[[1]]
			tt2[names(ps.res2$ps.individual)] <- ps.res2$ps.individual
			ps.res2$ps.individual <- tt2
		}
	}
	
	BIC.res <- NULL
	if(do.BIC) { BIC.res <- mclustBIC(data=dd[ , c("ESR1", "ERBB2"), drop=FALSE], modelNames=c(model.name), G=1:10)[ ,model.name] }

	## subtypes
	sbt[names(class.ts)] <- sbtn[class.ts]
	sbt.proba[dimnames(emclust.ts$z)[[1]], ] <- emclust.ts$z
	## discriminate between luminal A and B using AURKA
	sbt2 <- sbt
	sbt2[sbt == sbtn[3]] <- NA
	sbt2[sbt == sbtn[3] & dd2[ , "AURKA"] >= cc & complete.cases(sbt, dd2[ , "AURKA"])] <- sbtn2[3]
	sbt2[sbt == sbtn[3] & dd2[ , "AURKA"] < cc & complete.cases(sbt, dd2[ , "AURKA"])] <- sbtn2[4]	
	
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
		gplots::smartlegend(x="left", y="top", col=c("darkred", "darkgreen", "darkorange", "darkviolet"), legend=sbtn2, pch=c(17, 0, 10, 10), bg="white")
	}

	return(list("subtype"=sbt, "subtype.proba"=sbt.proba, "prediction.strength"=ps.res, "BIC"=BIC.res, "subtype2"=sbt2, "prediction.strength2"=ps.res2, "module.scores"=dd2, "mapping"=mymap))
}