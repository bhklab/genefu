#' @title Function for sample size calculation for correlation coefficients
#'
#' @description
#' This function enables to compute the sample size requirements for estimating 
#'   pearson, kendall and spearman correlations
#'
#' @usage
#' power.cor(rho, w, alpha = 0.05, method = c("pearson", "kendall", "spearman"))
#'
#' @param rho	Correaltion coefficients rho (Pearson, Kendall or Spearman)
#' @param w	a numerical vector of weights of the same length as x giving the weights to 
#'   use for elements of x in the first class.
#' @param alpha	alpha level
#' @param method	a character string specifying the method to compute the correlation 
#'   coefficient, must be one of "pearson" (default), "kendall" or "spearman". You can 
#'   specify just the initial letter.
#'
#' @return
#' sample size requirement
#'
#' @references
#' Bonett, D. G., and Wright, T. A. (2000). Sample size requirements for estimating 
#'   pearson, kendall and spearman correlations. Psychometrika, 65(1), 
#'   23-28. doi:10.1007/BF02294183
#'
#' @examples
#' power.cor(rho=0.5, w=0.1, alpha=0.05, method="spearman")
#'
#' @md
#' @export
#' 
## sample size calculation for correlation coefficients (Pearson, kendall and SPearman)
## example: power.cor(rho=0.5, w=0.1, alpha=0.05, method="spearman")
power.cor <- 
function (rho, w, alpha=0.05, method=c("pearson", "kendall", "spearman")) {
  method <- match.arg(method)
  bb <- c(3, 4, 3)
  cc <- c(1, sqrt(0.437), sqrt(1 + (rho^2 / 2)))
  names(bb) <- names(cc) <- c("pearson", "kendall", "spearman")
  bb <- bb[method]
  cc <- cc[method]
  nn0 <- 4 * cc^2 * (1 - rho^2)^2 * (qnorm(p=alpha/2, lower.tail=FALSE) / w) + bb
  if(nn0 < 10) { nn0t <- 10 } else { nn0t <- ceiling(nn0) }
  w0w <- sqrt(nn0t - bb) / sqrt(nn0 - bb)
  nn <- ceiling((nn0 - bb) * w0w^2 + bb)
  return(nn)
}