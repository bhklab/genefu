## weighted t test (Welch's version: unequal sample size, unequal variance)
## sources:
## http://en.wikipedia.org/wiki/T_test
## http://www.nicebread.de/blog/files/fc02e1635792cb0f2b3cbd1f7e6c580b-10.php
`weighted.ttest` <- 
function(x, w1, w2, alternative=c("two.sided", "less", "greater"), na.rm=FALSE) {
	alternative <- match.arg(alternative)
	ii <- complete.cases(x, w1, w2)
	if(!na.rm && sum(!ii) > 0) { stop("missing values are present!") } else { 
		w1 <- w1[ii]
		w2 <- w2[ii]
		x <- x[ii] 
	}
	if(!all(w1 >= 0 & w1 <= 1) || !all(w2 >= 0 & w2 <= 1) || (!all((w1 + w2) >= 0) && !all((w1 + w2) <= 1))) { stop("weights and their sum should lay in [0, 1]!") }
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