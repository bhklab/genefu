readarray<-function(dataFile,designFile=NA,hr=1,impute=T,method="mean"){
  
  headerRows <- hr
  
  x<-read.table(dataFile,sep="\t",header=F,fill=T,stringsAsFactors=FALSE)
  
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
    x<-read.table(designFile,sep="\t",header=T,row.names=1,fill=T,,stringsAsFactors=FALSE)
    xd<-xd[,sort.list(colnames(xd))]
    xd<-xd[,colnames(xd) %in% rownames(x)]
    x<-x[rownames(x) %in% colnames(xd),]
    x<-x[sort.list(rownames(x)),]
    classes<-as.data.frame(x)
  }
  
  if(sum(apply(xd,2,is.na))>0 & impute){
    library(impute)
    allAnn<-dimnames(xd)
    data.imputed<-impute.knn(as.matrix(xd))$data
    xd<-data.imputed[1:features,]
    dimnames(xd)<-allAnn
  }
  
  return(list(xd=xd, classes=classes, nfeatures=features, nsamples=samples, fnames=geneNames, snames=sampleNames))
}
