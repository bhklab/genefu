`npi` <-
function(size, grade, node, na.rm=FALSE) {

	nn <- names(size)
	if(is.null(nn)) { nn <- paste("X", 1:length(size), sep=".") }
	cc.ix <- complete.cases(size, grade, node)
	if(all(!cc.ix)) {
		tt <- rep(NA, length(size))
		names(tt) <- nn
		return(list("score"=tt, "risk"=tt))
	}
	size <- size[cc.ix]
	grade <- grade[cc.ix]
	node <- node[cc.ix]
	
	if(length(size) != length(grade) || length(grade) != length(node)) {
		stop("size, grade and lymph node stage must have the same length!")
	}
	if(!all(cc.ix) & !na.rm)  { stop("NA values are present!") }
	if(!all(is.element(grade, c("1", "2", "3")))) {
		stop("grade must be 1, 2 or 3!")
	}
	if(!all(is.element(node, c("1", "2", "3")))) {
		#if only "0" and "1" are available, map "0" -> "1" and "1" -> "3"
		stop("lymph node stage must be 1, 2 or 3!")
	}
	if(!is.numeric(size)) {
		stop("tumor size (cm) must be numeric!")
	}
	
	npi <- 0.2 * size + grade + node
	names(npi) <- nn[cc.ix]
	
	npi.score <- rep(NA, length(cc.ix))
	names(npi.score) <- nn
	npi.score[names(npi)] <- npi
	
	npi.c <- npi
	npi.c[npi < 3.4] <- "Good"
	npi.c[npi > 5.4] <- "Poor"
	npi.c[npi >= 3.4 & npi <= 5.4] <- "Intermediate"
	
	npi.classif <- rep(NA, length(cc.ix))
	names(npi.classif) <- nn
	npi.classif[names(npi.c)] <- npi.c
	#npi.classif <- as.factor(npi.classif)
	
	return(list("score"=npi.score, "risk"=npi.classif))
}