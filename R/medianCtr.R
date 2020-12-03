#' @title Center around the median
#'
#' @description
#' Utility function called within the claudinLow classifier
#'
#' @usage
#' medianCtr(x)
#'
#' @param x	 Matrix of numbers
#'
#' @return
#' A matrix of median-centered numbers
#' 
#' @references
#' citation("claudinLow")
#'
#' @seealso
#' [genefu::claudinLow]
#'
#' @md
#' @export
medianCtr<-function(x){
  annAll <- dimnames(x)
  medians <- apply(x,1,median,na.rm=T)
  x <- t(scale(t(x),center=medians,scale=F))
  dimnames(x) <- annAll
  return(x)
}
