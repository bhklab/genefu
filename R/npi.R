#' @title Function to compute the Nottingham Prognostic Index
#'
#' @description
#' This function computes the Nottingham Prognostic Index (NPI) as published
#'   in Galeat et al, 1992. NPI is a clinical index shown to be highly prognostic 
#'   in breast cancer.
#'
#' @usage
#' npi(size, grade, node, na.rm = FALSE)
#'
#' @param size tumor size in cm.
#' @param grade	Histological grade, i.e. low (1), intermediate (2) and high (3) grade.
#' @param node Nodal status. If only binary nodal status (0/1) is available, 
#'   map 0 to 1 and 1 to 3.
#' @param na.rm	TRUE if missing values should be removed, FALSE otherwise.
#'
#' @details
#' The risk prediction is either Good if score < 3.4, Intermediate 
#'   if 3.4 <= score <- 5.4, or Poor if score > 5.4.
#'
#' @return
#' A list with items:
#' - score: Continuous signature scores
#' - risk: Binary risk classification, 1 being high risk and 0 being low risk.


#' @references
#' Galea MH, Blamey RW, Elston CE, and Ellis IO (1992) "The nottingham 
#'   prognostic index in primary breast cancer", Breast Cancer Reasearch 
#'   and Treatment, 22(3):207-219.
#'
#' @seealso
#' [genefu::st.gallen]
#'
#' @examples
#' # load NKI dataset
#' data(nkis)
#' # compute NPI score and risk classification
#' npi(size=demo.nkis[ ,"size"], grade=demo.nkis[ ,"grade"],
#'   node=ifelse(demo.nkis[ ,"node"] == 0, 1, 3), na.rm=TRUE)
#' 
#' @md
#' @export
npi <-
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