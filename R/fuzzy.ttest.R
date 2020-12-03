#' @title Function to compute the fuzzy Student t test based on weighted 
#'   mean and weighted variance
#'
#' @description
#' This function allows for computing the weighted mean and weighted variance
#'   of a vector of continuous values.
#'
#' @usage
#' fuzzy.ttest(x, w1, w2, alternative=c("two.sided", "less", "greater"), 
#'   check.w = TRUE, na.rm = FALSE)
#'
#' @param x an object containing the observed values.
#' @param w1 a numerical vector of weights of the same length as x giving the weights
#'   to use for elements of x in the first class.
#' @param w2 a numerical vector of weights of the same length as x giving the weights to
#'   use for elements of x in the second class.
#' @param alternative a character string specifying the alternative hypothesis, must be one
#'   of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
#' @param check.w TRUE if weights should be checked such that 0 <= w <= 1 and (w1[i] + w2[i]) < 1 
#'   for 1 <= i <= length(x), FALSE otherwise. Beware that weights greater than one 
#'   may inflate over-optimistically resulting p-values, use with caution.
#' @param na.rm TRUE if missing values should be removed, FALSE otherwise.
#'
#' @details
#' The weights w1 and w2 should represent the likelihood for each observation stored in
#'   x to belong to the first and second class, respectively. Therefore the values contained
#'   in w1 and w2 should lay in [0,1] and 0 <= (w1[i] + w2[i]) <= 1 for i in {0,1,...,n} where
#'   n is the length of x.
#' The Welch's version of the t test is implemented in this function, therefore assuming
#'   unequal sample size and unequal variance. The sample size of the first and second class
#'   are calculated as the sum(w1) and sum(w2), respectively.
#'
#' @return
#' A numeric vector of six values that are the difference between the two weighted means, 
#'   the value of the t statistic, the sample size of class 1, the sample size of class 2, 
#'   the degree of freedom and the corresponding p-value.
#'
#' @references
#' http://en.wikipedia.org/wiki/T_test
#'
#' @seealso
#' [stats::weighted.mean]
#'
#'@examples
#' set.seed(54321)
#' # random generation of 50 normally distributed values for each of the two classes
#' xx <- c(rnorm(50), rnorm(50)+1)
#' # fuzzy membership to class 1
#' ww1 <- runif(50) + 0.3
#' ww1[ww1 > 1] <- 1
#' ww1 <- c(ww1, 1 - ww1)
#' # fuzzy membership to class 2
#' ww2 <- 1 - ww1
#' # Welch's t test weighted by fuzzy membership to class 1 and 2 
#' wt <- fuzzy.ttest(x=xx, w1=ww1, w2=ww2)
#' print(wt)
#' # Not run: 
#' # permutation test to compute the null distribution of the weighted t statistic
#' wt <- wt[2]
#' rands <- t(sapply(1:1000, function(x,y) { return(sample(1:y)) }, y=length(xx)))
#' randst <- apply(rands, 1, function(x, xx, ww1, ww2) 
#' { return(fuzzy.ttest(x=xx, w1=ww1[x], w2=ww2[x])[2]) }, xx=xx, ww1=ww1, ww2=ww2)
#' ifelse(wt < 0, sum(randst <= wt), sum(randst >= wt)) / length(randst) 
#' # End(Not run)
#'
#' @md
#' @export
fuzzy.ttest <- 
function(x, w1, w2, alternative=c("two.sided", "less", "greater"), check.w=TRUE, na.rm=FALSE) {
	alternative <- match.arg(alternative)
	ii <- complete.cases(x, w1, w2)
	if(!na.rm && sum(!ii) > 0) { stop("missing values are present!") } else { 
		w1 <- w1[ii]
		w2 <- w2[ii]
		x <- x[ii] 
	}
	if(check.w && (!all(w1 >= 0 & w1 <= 1) || !all(w2 >= 0 & w2 <= 1) || (!all((w1 + w2) >= 0) && !all((w1 + w2) <= 1)))) { stop("weights and their sum should lay in [0, 1]!") }
	tt <- weighted.meanvar(x=x, w=w1, na.rm=na.rm)
	x1.w <- tt[1]
	var1.w <- tt[2]
	tt <- weighted.meanvar(x=x, w=w2, na.rm=na.rm)
	x2.w <- tt[1]
	var2.w <- tt[2]
	n1 <- sum(w1)
	n2 <- sum(w2)
	t.value <- (x1.w - x2.w) / sqrt((var1.w / n1) + (var2.w / n2)) 
	df <- (((var1.w / n1) + (var2.w / n2))^2) / ((((var1.w / n1)^2) / (n1 - 1)) + (((var2.w / n2)^2) / (n2 - 1)))
	p.value <- pt(q=abs(t.value), df=df, lower.tail=FALSE)
	if(alternative == "two.sided")  { p.value <- p.value*2 }
	if(alternative == "less") { p.value <- 1-p.value }
	res <- c(x1.w - x2.w, t.value, n1, n2, df, p.value)
	names(res) <- c("diff", "t.value", "n1", "n2", "df", "p.value")
	return(res)
}