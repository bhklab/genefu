`strescR` <-
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