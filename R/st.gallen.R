#' @title Function to compute the St Gallen consensus criterion for 
#'   prognostication
#'
#' @description
#' This function computes the updated St Gallen consensus criterions as 
#'   published by Goldhirsh et al 2003.
#'
#' @usage
#' st.gallen(size, grade, node, her2.neu, age, vascular.inv, na.rm = FALSE)
#' 
#' @param size tumor size in cm.
#' @param grade Histological grade, i.e. low (1), intermediate (2) and 
#'   high (3) grade.
#' @param node Nodal status (0 or 1 for no lymph node invasion a,d at 
#'   least 1 invaded lymph ode respectively).
#' @param her2.neu Her2/neu status (0 or 1).
#' @param age Age at diagnosis (in years).
#' @param vascular.inv Peritumoral vascular invasion (0 or 1).
#' @param na.rm	TRUE if missing values should be removed, FALSE otherwise.
#'
#' @return
#' Vector of risk predictions: "Good", "Intermediate", and "Poor".
#'
#' @references
#' Goldhirsh A, Wood WC, Gelber RD, Coates AS, Thurlimann B, and Senn HJ 
#'   (2003) "Meeting highlights: Updated international expert 
#'   consensus on the primary therapy of early breast cancer", Journal of 
#'   Clinical Oncology, 21(17):3357-3365.
#'
#' @seealso
#' [genefu::npi]
#' 
#' @examples
#' # load NKI dataset
#' data(NKI)
#' # compute St Gallen predictions
#' st.gallen(size=demo.nkis[ ,"size"], grade=demo.nkis[ ,"grade"],
#'   node=demo.nkis[ ,"node"], her2.neu=sample(x=0:1, size=nrow(demo.nkis),
#'   replace=TRUE), age=demo.nkis[ ,"age"], vascular.inv=sample(x=0:1,                                                                                                           
#'   size=nrow(demo.nkis), replace=TRUE), na.rm=TRUE)
#'
#' @md
#' @export
st.gallen <-
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