`genius` <-
function(data, annot, do.mapping=FALSE, mapping, do.scale=TRUE) {
	
	## predict breast cancer molecular subtypes
	sbt.id <- subtype.cluster.predict(sbt.model=scmod1, data=data, annot=annot, do.mapping=do.mapping, mapping=mapping, do.scale=do.scale, do.prediction.strength=FALSE, do.BIC=FALSE, plot=FALSE, verbose=FALSE)
	
	usbt <- unique(sbt.id$subtype)
	usbt <- sort(usbt[!is.na(usbt)])
	pred.sbtclassif <- NULL
	for(ii in 1:length(usbt)) {
		myx <- sbt.id$subtype == usbt[ii] & !is.na(sbt.id$subtype)
		
		#compute the score from model
		score <- sig.score(x=sig.genius[[ii]][ , c("probe", "EntrezGene.ID",  "coefficient")], data=data, annot=annot, do.mapping=do.mapping,  mapping=mapping, verbose=FALSE)$score
		if(do.scale) {
			#the rescaling needs a large sample size!
			#necessary if we want to validate the classifier using a different dataset
			#the estimation of survival probabilities depends on the scale of the score
			score <-  (rescale(score, q=0.05, na.rm=TRUE) - 0.5) * 2
		}
		names(score) <- dimnames(data)[[1]]
		pred.sbtclassif <- c(pred.sbtclassif, list("score"=score))
	}
	names(pred.sbtclassif) <- names(sig.genius)
	#combine classifications
	cc <- NULL
	for(j in 1:length(pred.sbtclassif)) {
		cc <- cbind(cc, pred.sbtclassif[[j]])
	}
	ww <- sbt.id$subtype.proba
	combine.pred <- apply(ww * cc, 1, sum)

	pred.sbtclassif <- c(pred.sbtclassif, list(combine.pred))
	names(pred.sbtclassif)[length(pred.sbtclassif)] <- "score"
	return(pred.sbtclassif)
}