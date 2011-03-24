`read.m.file` <-
function(file, ...) {
	obj.file <- read.csv(file, stringsAsFactors=FALSE, ...)
	if(sum(!is.na(obj.file[ ,1]) & obj.file[ ,1] != "") > 0) {
		ix.delim <- c(which(obj.file[ ,1] != "")[-1]-1, nrow(obj.file) + 1)
		ix.f <- ix.l <- 1
		groups <- NULL
		npp <- np <- NULL
		for (i in 1:length(ix.delim)) {
			ix.l <- ix.delim[i] - 1
			np <- c(np, as.character(obj.file[ix.f,1]))
			groups <- c(groups, rep(i, ix.l - ix.f + 1))
			npp <- rbind(npp, obj.file[ix.f:ix.l,2:ncol(obj.file)])
			ix.f <- ix.l + 2
		}
		ugroups <- unique(groups)
		obj <- NULL
		for (j in 1:length(ugroups)) {
			obj <- c(obj, list(npp[groups == ugroups[j], ]))
		}
	names(obj) <- np
	} else { obj <- list("module"=obj)}
	return(obj)
}