`cordiff.dep` <-
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