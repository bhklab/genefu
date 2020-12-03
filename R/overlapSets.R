#' @title Overlap two datasets
#'
#' @description
#' Utility function called within the claudinLow classifien.
#'
#' @usage
#' overlapSets(x,y)
#'
#' @param x	Matrix1
#' @param y	Matrix2
#'
#' @return
#' A list of overlapped dataset
#'
#' @references
#' citation("claudinLow")
#'
#' @seealso
#' [genefu::claudinLow]
#'
#' @md
#' @export
overlapSets<-function(x,y){
  
  # subset the two lists to have a commonly ordered gene list
  x<-x[dimnames(x)[[1]] %in% dimnames(y)[[1]],]
  y<-y[dimnames(y)[[1]] %in% dimnames(x)[[1]],]
  
  #and sort such that thing are in the correct order
  x<-x[sort.list(row.names(x)),]
  y<-y[sort.list(row.names(y)),]
  
  return(list(x=x,y=y))
}
