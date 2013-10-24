`spearmanCI` <- 
function (x, n, alpha=0.05) {
    zz <- sqrt((n-3)/1.06) * survcomp::fisherz(x)
    zz.se <- 1/sqrt(n - 3)
    ll <- zz - qnorm(p=alpha/2, lower.tail=FALSE) * zz.se
    ll <- survcomp::fisherz(ll / sqrt((n-3)/1.06), inv=TRUE)
    uu <- zz + qnorm(p=alpha/2, lower.tail=FALSE) * zz.se
    uu <- survcomp::fisherz(uu/ sqrt((n-3)/1.06), inv=TRUE)
    pp <- pnorm(q=zz, lower.tail=x<0)
    res <- c("lower"=ll, "upper"=uu, "p.value"=pp)
    return(res)
}