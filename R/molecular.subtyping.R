if(getRversion() >= "2.15.1")  utils::globalVariables(c("scmgene.robust","scmod2.robust","pam50.robust","ssp2006.robust","ssp2003.robust","claudinLowData"))

molecular.subtyping <- 
function (sbt.model=c("scmgene", "scmod1", "scmod2", "pam50", "ssp2006", "ssp2003", "intClust", "AIMS","claudinLow"), data, annot, do.mapping=FALSE) {
  
  sbt.model <- match.arg(sbt.model)
  
  ## convert SCM to SSP nomenclature
  sbt.conv <- rbind(c("ER-/HER2-", "Basal"),
    c("HER2+", "Her2"),
    c("ER+/HER2- High Prolif", "LumB"),
    c("ER+/HER2- Low Prolif", "LumA")
  )
  colnames(sbt.conv) <- c("SCM.nomenclature", "SSP.nomenclature")
    
  sbtn.ssp <- c("Basal", "Her2", "LumB", "LumA", "Normal")
  sbtn2.ssp <- c("Basal", "Her2", "Lums", "LumB", "LumA", "Normal")
  
  ## SCM family
  if (sbt.model %in% c("scmgene", "scmod1", "scmod2")) {
    switch(sbt.model,
      "scmgene" = {
        sbts <- subtype.cluster.predict(sbt.model=scmgene.robust, data=data, annot=annot, do.mapping=do.mapping)[c("subtype2", "subtype.proba2")]
      },
      "scmod1" = {
        sbts <- subtype.cluster.predict(sbt.model=scmod1.robust, data=data, annot=annot, do.mapping=do.mapping)[c("subtype2", "subtype.proba2")]
      },
      "scmod2" = {
        sbts <- subtype.cluster.predict(sbt.model=scmod2.robust, data=data, annot=annot, do.mapping=do.mapping)[c("subtype2", "subtype.proba2")]
      }
    )
    names(sbts) <- c("subtype", "subtype.proba")
    ## compute crisp classification
    sbts$subtype.crisp <- t(apply(sbts$subtype.proba, 1, function (x) {
      xx <- array(0, dim=length(x), dimnames=list(names(x)))
      xx[which.max(x)] <- 1
      return (xx)
    }))
    
    ## reorder columns
    #ss <- sbtn2.ssp[is.element(sbtn2.ssp, colnames(sbts$subtype.proba))]
    #sbts$subtype.proba <- sbts$subtype.proba[ , ss, drop=FALSE]
    #sbts$subtype.crisp <- sbts$subtype.crisp[ , ss, drop=FALSE]
    
    ## set the proper names
    names(sbts$subtype) <- rownames(sbts$subtype.proba) <- rownames(sbts$subtype.crisp)<- rownames(data)
  }
  
  ## SSP family
  if (sbt.model %in% c("ssp2003", "ssp2006", "pam50")) {
    switch(sbt.model,
      "pam50" = {
        sbts <- intrinsic.cluster.predict(sbt.model=pam50.robust, data=data, annot=annot, do.mapping=do.mapping)[c("subtype", "subtype.proba")]
      },
      "ssp2006" = {
        sbts <- intrinsic.cluster.predict(sbt.model=ssp2006.robust, data=data, annot=annot, do.mapping=do.mapping)[c("subtype", "subtype.proba")]
      },
      "ssp2003" = {
        sbts <- intrinsic.cluster.predict(sbt.model=ssp2003.robust, data=data, annot=annot, do.mapping=do.mapping)[c("subtype", "subtype.proba")]
      }
    )
    sbts$subtype <- factor(as.character(sbts$subtype), levels=sbtn.ssp)
    ## compute crisp classification
    sbts$subtype.crisp <- t(apply(sbts$subtype.proba, 1, function (x) {
      xx <- array(0, dim=length(x), dimnames=list(names(x)))
      xx[which.max(x)] <- 1
      return (xx)
    }))
    
    ## merge LumA and LumB: #sum the probability for LumA and LumB to get the probability for Luminals in general
    #lums.proba <- apply(sbts$subtype.proba[ , c("LumB", "LumA"), drop=FALSE], 1, sum, na.rm=TRUE)
    #sbts$subtype.proba <- cbind(sbts$subtype.proba, "Lums"=lums.proba)
    #lums.crisp <- as.numeric(is.element(sbts$subtype, c("LumA", "LumB")))
    #sbts$subtype.crisp <- cbind(sbts$subtype.crisp, "Lums"=lums.crisp)
    
    ## reorder columns
    #ss <- sbtn2.ssp[is.element(sbtn2.ssp, colnames(sbts$subtype.proba))]
    #sbts$subtype.proba <- sbts$subtype.proba[ , ss, drop=FALSE]
    #sbts$subtype.crisp <- sbts$subtype.crisp[ , ss, drop=FALSE]
    
    ## set the proper names
    names(sbts$subtype) <- rownames(sbts$subtype.proba) <- rownames(sbts$subtype.crisp)<- rownames(data)
  }
  
  ## IntClust family
  if (sbt.model %in% c("intClust")) {
    #message("Note: Need a Gene.Symbol column in the annotation object")
    sbts<-NULL
    myx <- !is.na(annot[ , "Gene.Symbol"]) & !duplicated(annot[ , "Gene.Symbol"])
    dd <- t(data[ , myx, drop=FALSE])
    rownames(dd) <- annot[myx, "Gene.Symbol"]
    ## remove patients with more than 80% missing values
    rix <- apply(dd, 2, function (x, y) { return ((sum(is.na(x) / length(x))) > y) }, y=0.8)
    cix <- apply(dd, 2, function (x, y) { return ((sum(is.na(x) / length(x))) > y) }, y=0.8)
    dd <- dd[!rix, !cix, drop=FALSE]
    features <- iC10::matchFeatures(Exp=dd, Exp.by.feat="Gene.Symbol")
    features <- iC10::normalizeFeatures(features, method="scale")
    res <- iC10::iC10(features)
    ## compute crisp classification
    crisp <- t(apply(res$posterior, 1, function (x) {
      xx <- array(0, dim=length(x), dimnames=list(names(x)))
      xx[which.max(x)] <- 1
      return (xx)
    }))
    sbts$subtype <- array(NA, dim=nrow(data), dimnames=list(rownames(data)))
    sbts$subtype[!rix] <- res$class
    sbts$subtype.proba <- array(NA, dim=c(nrow(data), ncol(res$posterior)), dimnames=list(rownames(data), colnames(res$posterior)))
    sbts$subtype.proba[!rix, ] <- res$posterior
    sbts$subtype.crisp <- t(apply(sbts$subtype.proba, 1, function (x) {
      xx <- array(0, dim=length(x), dimnames=list(names(x)))
      xx[which.max(x)] <- 1
      return (xx)
    }))
    ## set the proper colnames
    colnames(sbts$subtype.proba) <- colnames(sbts$subtype.crisp) <- paste("iC", colnames(sbts$subtype.proba), sep="")
    sbts$subtype <- paste("iC", sbts$subtype, sep="")
    sbts$subtype <- factor(sbts$subtype, levels=colnames(sbts$subtype.proba))
    ## set the proper rownames
    names(sbts$subtype) <- rownames(sbts$subtype.proba) <- rownames(sbts$subtype.crisp)<- rownames(data)
  }
  
  ## AIMS classifier
  if (sbt.model %in% c("AIMS")) {
      sbts <- AIMS::applyAIMS(eset=t(data), EntrezID=annot[ , "EntrezGene.ID"])[c("cl", "all.probs")]
      sbts$subtype <- sbts$cl
      sbts$subtype.proba <- matrix(unlist(sbts$all.probs$`20`), ncol = 5, byrow = TRUE)
      colnames(sbts$subtype.proba) <- colnames(sbts$all.probs$`20`)
      rownames(sbts$subtype.proba) <- rownames(sbts$subtype)

      ## compute crisp classification
      sbts$subtype.crisp <- t(
        apply(sbts$subtype.proba, 1, function (x) {
        xx <- array(0, dim=length(x), dimnames=list(names(x)))
        xx[which.max(x)] <- 1
        return (xx)
      })
      )
      sbts<-sbts[- which(names(sbts) %in% c("cl","all.probs"))]
  }
  
  ## CLAUDIN-LOW classifier
  if (sbt.model %in% c("claudinLow")) {
    train<-claudinLowData
    train$xd<- medianCtr(train$xd)
    
    if(do.mapping) {
      gid1 <- as.numeric(rownames(train$xd))
      names(gid1) <- rownames(train$xd)
      gid2 <- as.numeric(as.character(annot[ ,"EntrezGene.ID"]))
      names(gid2) <- dimnames(annot)[[1]]

      ## remove missing and duplicated geneids from the gene list
      rm.ix <- is.na(gid1) | duplicated(gid1)
      gid1 <- gid1[!rm.ix]
      rr <- geneid.map(geneid1=gid2, data1=data, geneid2=gid1, verbose=FALSE)
      gm <- length(rr$geneid2)
      if(is.na(rr$geneid1[1])) {
        gm <- 0
        #no gene ids in common
        res <- rep(NA, nrow(data))
        names(res) <- dimnames(data)[[1]]
        gf <- c("mapped"=0, "total"=gt)
        if(verbose) { message(sprintf("probe candidates: 0/%i", gt)) }
        return(list("score"=res, "risk"=res, "mapping"=gf, "probe"=NA))
      }
      gid1 <- rr$geneid2
      gid2 <- rr$geneid1
      data <- rr$data1
      #mymapping <- c("mapped"=gm, "total"=gt)
      myprobe <- cbind("probe"=names(gid1), "EntrezGene.ID"=gid1, "new.probe"=names(gid2))
      ## change the names of probes in the data
      dimnames(data)[[2]] <- names(gid2) <- names(gid1)
    } 
    
    test<- medianCtr(t(data)) #probes as rows, median-centered
    #Run Classifier Call
    sbts<-claudinLow(train$xd, as.matrix(train$classes$Group,ncol=1), test)
    sbts$subtype<-factor(sbts$predictions$Call)
    
#     sbts$subtype.proba <- t(apply(X=predout$centroids[,2], MARGIN=1, FUN=function(x) { return(abs(x) / sum(abs(x), na.rm=TRUE)) }))
#     
#     ## compute crisp classification
#     sbts$subtype.crisp <- t(
#       apply(sbts$subtype.proba, 1, function (x) {
#         xx <- array(0, dim=length(x), dimnames=list(names(x)))
#         xx[which.max(x)] <- 1
#         return (xx)
#       })
#     )
#     
  }
  return (sbts)
}