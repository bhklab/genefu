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
#' @param w	The desired confidence interval width
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
    za2 <- qnorm(p=alpha/2, lower.tail=FALSE)
    nn0 <- 4 * cc^2 * (1 - rho^2)^2 * (za2 / w)^2 + bb
    if(nn0 < 10) { nn0t <- 10 } else { nn0t <- ceiling(nn0) }

    L1 = 0.5*(log(1+rho) - log(1-rho)) - cc*za2/(sqrt(nn0t - bb))
    L2 = 0.5*(log(1+rho) - log(1-rho)) + cc*za2/(sqrt(nn0t - bb))
    
    lowerLimit <- (exp(2*L1)-1)/(exp(2*L1)+1)
    upperLimit <- (exp(2*L2)-1)/(exp(2*L2)+1)
    
    w0 <- abs(upperLimit-lowerLimit)
    
    w0w <- w0 / w
    nn <- ceiling((nn0 - bb) * w0w^2 + bb)
    
    return(nn)
  }
