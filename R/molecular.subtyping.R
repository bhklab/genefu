molecular.subtyping <- 
function (sbt.model=c("scmgene", "scmod1", "scmod2", "pam50", "ssp2006", "ssp2003", "intClust", "AIMS"), data, annot, do.mapping=FALSE) {
  
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
    ss <- sbtn2.ssp[is.element(sbtn2.ssp, colnames(sbts$subtype.proba))]
    sbts$subtype.proba <- sbts$subtype.proba[ , ss, drop=FALSE]
    sbts$subtype.crisp <- sbts$subtype.crisp[ , ss, drop=FALSE]
    
    ## set the proper names
    names(sbts$subtype) <- rownames(sbts$subtype.proba) <- rownames(sbts$subtype.crisp)<- rownames(data)
  }
  
  ## IntClust family
  if (sbt.model %in% c("intClust")) {
    
    myx <- !is.na(annot[ , "Gene.Symbol"]) & !duplicated(annot[ , "Gene.Symbol"])
    dd <- t(data[ , myx, drop=FALSE])
    rownames(dd) <- annot[myx, "gene"]
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
      AIMS::applyAIMS(eset=t(data), EntrezID=annot[ , "EntrezGene.ID"])
  }
  
  return (sbts)
}