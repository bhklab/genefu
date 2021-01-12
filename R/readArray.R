#' @title Overlap two datasets
#'
#' @description
#' Formatting function to read arrays and format for use in the claudinLow classifier.
#'
#' @usage
#' readArray(dataFile,designFile=NA,hr=1,impute=TRUE,method="mean")
#'
#' @param dataFile file with matrix to be read.
#' @param designFile Design of file.
#' @param hr Header rows as Present (2) or Absent (1).
#' @param impute whether data will be imputed or not.
#' @param method Default method is "mean".
#'
#' @return
#' A list
#'
#' @references
#' citation("claudinLow")
#'
#' @seealso
#' [genefu::claudinLow]
#'
#' @md
#' @importFrom impute impute.knn
#' @export
readArray <- function(dataFile, designFile=NA, hr=1, impute=TRUE,
    method="mean")
{

  headerRows <- hr

  x<-read.table(dataFile,sep="\t",header=FALSE,fill=TRUE,stringsAsFactors=FALSE)

  if(headerRows==1){
    sampleNames<-as.vector(t(x[1,-1]))
    x<-x[-1,]
    classes<-NULL
    ids<-x[,1]
    xd<-x[,-1]
    xd<-apply(xd,2,as.numeric)
    xd<-collapseIDs(xd,ids,method)
  }else{
    sampleNames<-as.vector(t(x[1,-1]))
    x<-x[-1,]

    classes<-x[1:(headerRows-1),]
    dimnames(classes)[[1]]<-classes[,1]
    classes<-classes[,-1]
    classes[classes==""]<-NA
    classes<-t(classes)
    rownames(classes)<-sampleNames
    classes<-as.data.frame(classes)

    xd<-x[(-1:-(headerRows-1)),]
    ids<-as.vector(t(xd[,1]))
    xd<-xd[,-1]
    xd<-apply(xd,2,as.numeric)
    xd<-collapseIDs(xd,ids,method)
  }

  features<- dim(xd)[1]
  samples<- dim(xd)[2]
  geneNames<-rownames(xd)
  xd<-apply(xd,2,as.numeric)
  rownames(xd)<-geneNames
  colnames(xd)<-sampleNames

  if(!is.na(designFile)){
    x<-read.table(designFile,sep="\t", header=TRUE, row.names=1, fill=TRUE,
                  stringsAsFactors=FALSE)
    xd<-xd[,sort.list(colnames(xd))]
    xd<-xd[,colnames(xd) %in% rownames(x)]
    x<-x[rownames(x) %in% colnames(xd),]
    x<-x[sort.list(rownames(x)),]
    classes<-as.data.frame(x)
  }

  if(sum(apply(xd,2,is.na))>0 & impute){
    #library(impute)
    allAnn<-dimnames(xd)
    data.imputed<-impute.knn(as.matrix(xd))$data
    xd<-data.imputed[1:features,]
    dimnames(xd)<-allAnn
  }

  return(list(xd=xd, classes=classes, nfeatures=features, nsamples=samples, fnames=geneNames, snames=sampleNames))
}
