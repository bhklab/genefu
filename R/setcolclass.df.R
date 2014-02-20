`setcolclass.df` <-
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