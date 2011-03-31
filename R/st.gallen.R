`st.gallen` <-
function(size, grade, node, her2.neu, age, vascular.inv, na.rm=FALSE) {

	nn <- names(size)
	if(is.null(nn)) { nn <- paste("PATIENT",  1:length(size),  sep=".") }
	names(size) <- names(grade) <- names(node) <- names(her2.neu) <- names(age) <- names(vascular.inv) <- nn
	
	## remove missing values
	cc.ix <- complete.cases(size, grade, node, her2.neu, age, vascular.inv)
	size <- size[cc.ix]
	grade <- grade[cc.ix]
	node <- node[cc.ix]
	her2.neu <- her2.neu[cc.ix]
	age <- age[cc.ix]
	vascular.inv <- vascular.inv[cc.ix]
	
	if(length(size) + length(grade) + length(node) + length(her2.neu) + length(age) + length(vascular.inv) != (6 * length(size))) {
		stop("size, grade, lymph node stage, her2/neu expression, age and peritumoral vascular invasion must have the same length!")
	}
	if(!all(cc.ix) & !na.rm)  { stop("NA values are present!") }
	if(!all(is.element(grade, c("1", "2", "3")))) {
		stop("grade must be 1, 2 or 3!")
	}
	if(!all(is.element(node, c("0", "1")))) {
		stop("lymph node stage must be 0 or 1!")
	}
	if(!is.numeric(size)) {
		stop("tumor size (cm) must be numeric!")
	}
	if(!is.numeric(age)) {
		stop("age (years) must be numeric!")
	}
	if(!all(is.element(her2.neu, c("0", "1")))) {
		stop("her2/neu expression must be 0 or 1!")
	}
	if(!all(is.element(vascular.inv, c("0", "1")))) {
		stop("peritumoral vascular invasion must be 0 or 1!")
	}
	
	lowr <- node == 0 & (size <= 2 & grade == 1 & vascular.inv == 0 & her2.neu == 0 & age >= 35)
	names(lowr) <- nn[cc.ix]
	intermediater <- (node == 0 & (size > 2 | grade != 1 | vascular.inv == 1 | her2.neu == 1 | age < 35)) | (node == 1 & her2.neu == 0)
	names(intermediater) <- nn[cc.ix]
	highr <- (node == 1 & (her2.neu == 1)) # | (node.stage == 3)
	names(highr) <- nn[cc.ix]
	
	#if(sum(lowr) + sum(highr) + sum(intermediater) != (3 * sum(cc.ix))) {
	#	stop("the classification is not unique!")
	#}
	
	stgr <- rep(NA, length(cc.ix))
	names(stgr) <- nn
	stgr[names(lowr)][lowr] <- "Good"
	stgr[names(intermediater)][intermediater] <- "Intermediate"
	stgr[names(highr)][highr] <- "Poor"
	
	return(stgr)
}