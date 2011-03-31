`rename.duplicate` <-
function (x, sep="_", verbose=FALSE) {

	x <- as.character(x)
	duplix <- duplicated(x)
	duplin <- x[duplix]

	ix <- numeric(length=length(unique(duplin)))
	names(ix) <- unique(duplin)
	retval <- numeric(length=length(duplin))
	for(i in 1:length(duplin)) { retval[i] <- ix[duplin[i]] <- ix[duplin[i]] + 1 }
	retval <- retval + 1
	x[duplix] <- paste(duplin, retval, sep=sep)

	if (verbose) { message(sprintf("%i duplicated names", length(duplin))) }
	
	return (list(new.x=x, duplicated.x=duplin))
}