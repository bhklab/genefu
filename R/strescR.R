#' @title Utility function to escape LaTeX special characters 
#'   present in a string
#'
#' @description
#' This function returns a vector of strings in which LaTeX special 
#'   characters are escaped, this was useful in conjunction with xtable.
#'
#' @usage
#' strescR(strings)
#'
#' @param strings	A vector of strings to deal with.
#'
#' @return
#' A vector of strings with escaped characters within each string.
#'
#' @references
#' citation("seqinr")
#'
#' @seealso
#' stresc
#' 
#' @examples
#' strescR("MISC_RNA")
#' strescR(c("BB_0001","BB_0002"))
#'
#' @md
#' @export
strescR <-
function (strings) {

	c2s <- function (chars = c("m", "e", "r", "g", "e", "d")) {
		return(paste(chars, collapse = ""))
	}
	
	s2c <- function (string) {
		if (is.character(string) && length(string) == 1) {
			return(unlist(strsplit(string, split = "")))
		} else {
			warning("Wrong argument type in s2c(), NA returned")
			return(NA)
		}
	}

	fromchar <- s2c("\\{}$^_%#&~[]")
	tochar <- c("$\\backslash$", "\\{", "\\}", "\\$", "\\^{}", 
	"\\_", "\\%", "\\#", "\\&", "\\~{}", "\\[", "\\]")
	f <- function(string) {
		c2s(sapply(s2c(string), function(x) ifelse(x %in% fromchar, tochar[which(x == fromchar)], x)))
	}
	
	return(sapply(strings, f, USE.NAMES = FALSE))
}