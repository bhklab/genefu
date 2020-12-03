#' @title Function to set the class of columns in a data.frame
#'
#' @description
#' This function enables to set the class of each column in a data.frame.
#'
#' @usage
#' setcolclass.df(df, colclass, factor.levels)
#'
#' @param df data.frame for which columns' class need to be updated.
#' @param colclass class for each column of the data.frame.
#' @param factor.levels	list of levels for each factor.
#'
#' @return
#' A data.frame with columns' class and levels properly set
#'
#' @examples
#' tt <- data.frame(matrix(NA, nrow=3, ncol=3, dimnames=list(1:3, paste("column", 1:3, sep="."))), 
#'   stringsAsFactors=FALSE)
#' tt <- setcolclass.df(df=tt, colclass=c("numeric", "factor", "character"), 
#'   factor.levels=list(NULL, c("F1", "F2", "F3"), NULL))
#'
#' @md
#' @export
setcolclass.df <-
function (df, colclass, factor.levels) {
	ww <- options()$warn
	options(warn=-1)
	toCls <- function(x, cls) { do.call(paste("as", cls, sep = "."), list(x)) }
	df <- replace(df, , Map(toCls, x=df, cls=colclass))
	options(warn=ww)
	iix <- FALSE
	if(!missing(factor.levels)) { iix <- colclass == "factor" & !is.null(factor.levels) }
	if(any(iix)) {
		for(i in which(iix)) { levels(df[[i]]) <- factor.levels[[i]] }
	}
	return(df)
}