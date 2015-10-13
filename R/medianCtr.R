medianCtr<-function(x){
  annAll <- dimnames(x)
  medians <- apply(x,1,median,na.rm=T)
  x <- t(scale(t(x),center=medians,scale=F))
  dimnames(x) <- annAll
  return(x)
}
