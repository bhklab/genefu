#' @title Function to rename duplicated strings
#'
#' @description
#' This function renames duplicated strings by adding their number of 
#'   occurrences at the end.
#'
#' @usage
#' rename.duplicate(x, sep = "_", verbose = FALSE)
#'
#' @param x	vector of strings.
#' @param sep	a character to be the separator between the number added at 
#'   the end and the string itself.
#' @param verbose	TRUE to print informative messages, FALSE otherwise.
#'
#' @return
#' A list with items:
#' - new.x:	new strings (without duplicates).
#' - duplicated.x: strings which were originally duplicated.
#'
#' @examples
#' nn <- sample(letters[1:10], 30, replace=TRUE)
#' table(nn)
#' rename.duplicate(x=nn, verbose=TRUE)
#'
#' @md
#' @export
rename.duplicate <-
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