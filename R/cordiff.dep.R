#' @title Function to estimate whether two dependent correlations differ
#'
#' @description
#' This function tests for statistical differences between two dependent correlations
#'   using the formula provided on page 56 of Cohen & Cohen (1983). The function returns 
#'   a t-value, the DF and the p-value.
#'
#' @usage
#' cordiff.dep(r.x1y, r.x2y, r.x1x2, n,
#'   alternative = c("two.sided", "less", "greater"))
#'
#' @param r.x1y	The correlation between x1 and y where y is typically your outcome variable.
#' @param r.x2y	The correlation between x2 and y where y is typically your outcome variable.
#' @param r.x1x2 The correlation between x1 and x2 (the correlation between your two predictors).
#' @param n The sample size.
#' @param alternative A character string specifying the alternative hypothesis, must be
#'   one of "two.sided" default), "greater" or "less". You can specify just the initial letter.
#'
#' @details
#' This function is inspired from the cordif.dep.
#'
#' @return
#' Vector of three values: t statistics, degree of freedom, and p-value.
#'
#' @references
#' Cohen, J. & Cohen, P. (1983) "Applied multiple regression/correlation analysis for the
#'   behavioral sciences (2nd Ed.)" Hillsdale, nJ: Lawrence Erlbaum Associates.
#'
#' @seealso
#' [stats::cor], [stats::t.test], [genefu::compareProtoCor]
#'
#' @examples
#' # load VDX dataset
#' data(vdxs)
#' # retrieve ESR1, AURKA and MKI67 gene expressions
#' x1 <- data.vdxs[ ,"208079_s_at"]
#' x2 <- data.vdxs[ ,"205225_at"]
#' y <- data.vdxs[ ,"212022_s_at"]
#' # is MKI67 significantly more correlated to AURKA than ESR1?
#' cc.ix <- complete.cases(x1, x2, y)
#' cordiff.dep(r.x1y=abs(cor(x=x1[cc.ix], y=y[cc.ix], use="everything",
#'   method="pearson")), r.x2y=abs(cor(x=x2[cc.ix], y=y[cc.ix],
#'   use="everything", method="pearson")), r.x1x2=abs(cor(x=x1[cc.ix],
#'   y=x2[cc.ix], use="everything", method="pearson")), n=sum(cc.ix),
#'   alternative="greater")
#'
#' @md
#' @export
cordiff.dep <-
function(r.x1y, r.x2y, r.x1x2, n, alternative=c("two.sided", "less", "greater")) {
	alternative <- match.arg(alternative)
	rbar <- (r.x1y + r.x2y)/2
	barRbar <- 1 - r.x1y^2 - r.x2y^2 - r.x1x2^2 + 2 * r.x1y * r.x2y * r.x1x2
	tvalue.num <- ((r.x1y - r.x2y) * sqrt((n - 1) * (1 + r.x1x2)))
	tvalue.den <- sqrt(((2 * ((n - 1)/(n - 3))) * barRbar + ((rbar^2)) * (1 - r.x1x2)^3))
	t.value <- tvalue.num / tvalue.den
	DF <- n - 3
	switch(alternative,
	"greater"={
		p.value <- pt(t.value, DF, lower.tail=FALSE)
	},
	"less"={
		p.value <- 1 - pt(t.value, DF, lower.tail=FALSE)
	},
	"two.sided"={
		p.value <- (1 - pt(abs(t.value), DF)) * 2
	})
	OUT <- c(t.value, DF, p.value)
	names(OUT) <- c("t.value",  "DF",  "p.value")
	
	return(OUT)
}