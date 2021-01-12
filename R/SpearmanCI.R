#' @title Function to compute the confidence interval for the Spearman 
#'   correelation coefficient
#'
#' @description
#' This function enables to compute the confidence interval for the Spearman 
#'   correelation coefficient using the Fischer Z transformation.
#'
#' @usage
#' spearmanCI(x, n, alpha = 0.05)
#'
#' @param x	Spearman correlation coefficient rho.
#' @param n	the sample size used to compute the Spearman rho.
#' @param alpha	alpha level for confidence interval.
#'
#' @return
#' A vector containing the lower, upper values for the confidence interval 
#'   and p-value for Spearman rho
#'
#' @examples
#' spearmanCI(x=0.2, n=100, alpha=0.05)
#'
#' @md
#' @export
spearmanCI <- 
function (x, n, alpha=0.05) {
    zz <- sqrt((n-3)/1.06) * survcomp::fisherz(x)
    zz.se <- 1/sqrt(n - 3)
    ll <- zz - qnorm(p=alpha / 2, lower.tail=FALSE) * zz.se
    ll <- survcomp::fisherz(ll / sqrt((n-3)/1.06), inv=TRUE)
    uu <- zz + qnorm(p=alpha / 2, lower.tail=FALSE) * zz.se
    uu <- survcomp::fisherz(uu / sqrt((n-3)/1.06), inv=TRUE)
    pp <- pnorm(q=zz, lower.tail=x < 0)
    res <- c("lower"=ll, "upper"=uu, "p.value"=pp)
    return(res)
}