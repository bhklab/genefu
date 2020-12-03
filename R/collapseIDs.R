#' @title Utility function to collapse IDs
#'
#' @description
#' Utility function called within the claudinLow classifier
#'
#' @usage
#' collapseIDs(x,allids=row.names(x),method="mean")
#'
#' @param x Matrix of numbers.
#' @param allids Defaults to rownames of matrix.
#' @param method Default method is "mean".
#'
#'
#' @return
#' A matrix
#'
#' @references
#' citation("claudinLow")
#'
#' @seealso
#' [genefu::claudinLow]
#'
#' @md
#' @export
collapseIDs<-function(x,allids=row.names(x),method="mean"){
  
  allids<-as.vector(allids)
  ids<- levels(as.factor(allids))
  x.col<- matrix(nrow=length(ids), ncol=dim(x)[2])
  
  if(length(ids)==dim(x)[1]){ 
    dimnames(x)[[1]]<-allids
    return(x) 
  }
  
  for(i in 1:length(ids)){
    if(sum(allids==ids[i])>1){
      indices <- allids==ids[i] 
      if(method=="mean"){
        vals<-apply(x[indices,],2,mean,na.rm=T)
      }
      if(method=="median"){
        vals<-apply(x[indices,],2,median,na.rm=T)
      }
      if(method=="stdev"){   
        temp<- x[indices,]
        stdevs<- apply(temp,1,sd,na.rm=T)
        vals<- temp[match(max(stdevs),stdevs),]
      }
      if(method=="sum"){   
        vals<-apply(x[indices,],2,sum,na.rm=T)
      }
      if(method=="iqr"){   
        temp<- x[indices,]
        iqrs<- apply(temp,1,function(x){quantile(x,.75,na.rm=T)-quantile(x,.25,na.rm=T)})
        vals<- temp[match(max(iqrs),iqrs),]
      }
      x.col[i,] <- vals
    }else{
      x.col[i,] <- t(as.vector(x[allids==ids[i],]))
    }
  }
  
  dimnames(x.col)<- list(ids,dimnames(x)[[2]])
  return(x.col)
  
}
