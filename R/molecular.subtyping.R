if(getRversion() >= "2.15.1")
    utils::globalVariables(c("scmgene.robust","scmod2.robust","pam50.robust",
                             "ssp2006.robust","ssp2003.robust","claudinLowData"))

#' @title Function to identify breast cancer molecular subtypes using
#'   the Subtype Clustering Model
#'
#' @description
#' This function identifies the breast cancer molecular subtypes using a Subtype
#'   Clustering Model fitted by subtype.cluster.
#'
#' @usage
#' molecular.subtyping(sbt.model = c("scmgene", "scmod1", "scmod2",
#'   "pam50", "ssp2006", "ssp2003", "intClust", "AIMS","claudinLow"),
#'   data, annot, do.mapping = FALSE, verbose = FALSE)
#'
#' @param sbt.model	Subtyping classification model, can be either "scmgene", "scmod1",
#'   "scmod2", "pam50", "ssp2006", "ssp2003", "intClust", "AIMS", or "claudinLow".
#' @param data Matrix of gene expressions with samples in rows and probes in columns,
#'   dimnames being properly defined.
#' @param annot Matrix of annotations with at least one column named "EntrezGene.ID"
#'   (for ssp, scm, AIMS, and claudinLow models) or "Gene.Symbol" (for the intClust
#'   model), dimnames being properly defined.
#' @param do.mapping TRUE if the mapping through Entrez Gene ids must be performed
#'   (in case of ambiguities, the most variant probe is kept for each gene), FALSE otherwise.
#' @param verbose TRUE if informative messages should be displayed, FALSE otherwise.
#'
#' @return
#' A list with items:
#' - subtype: Subtypes identified by the subtyping classification model.
#' - subtype.proba: Probabilities to belong to each subtype estimated by the
#'   subtyping classification model.
#' - subtype.crisp: Crisp classes identified by the subtyping classification model.
#'
#' @references
#' T. Sorlie and R. Tibshirani and J. Parker and T. Hastie and J. S. Marron and A.
#'   Nobel and S. Deng and H. Johnsen and R. Pesich and S. Geister and J. Demeter and
#'   C. Perou and P. E. Lonning and P. O. Brown and A. L. Borresen-Dale and D. Botstein
#'   (2003) "Repeated Observation of Breast Tumor Subtypes in Independent Gene
#'   Expression Data Sets", Proceedings of the National Academy of Sciences,
#'   1(14):8418-8423
#' Hu, Zhiyuan and Fan, Cheng and Oh, Daniel and Marron, JS and He, Xiaping and
#'   Qaqish, Bahjat and Livasy, Chad and Carey, Lisa and Reynolds, Evangeline and
#'   Dressler, Lynn and Nobel, Andrew and Parker, Joel and Ewend, Matthew and Sawyer,
#'   Lynda and Wu, Junyuan and Liu, Yudong and Nanda, Rita and Tretiakova, Maria and
#'   Orrico, Alejandra and Dreher, Donna and Palazzo, Juan and Perreard, Laurent and
#'   Nelson, Edward and Mone, Mary and Hansen, Heidi and Mullins, Michael and
#'   Quackenbush, John and Ellis, Matthew and Olopade, Olufunmilayo and Bernard,
#'   Philip and Perou, Charles (2006) "The molecular portraits of breast tumors are
#'   conserved across microarray platforms", BMC Genomics, 7(96)
#' Parker, Joel S. and Mullins, Michael and Cheang, Maggie C.U. and Leung, Samuel and
#'   Voduc, David and Vickery, Tammi and Davies, Sherri and Fauron, Christiane and He,
#'   Xiaping and Hu, Zhiyuan and Quackenbush, John F. and Stijleman, Inge J. and Palazzo,
#'   Juan and Marron, J.S. and Nobel, Andrew B. and Mardis, Elaine and Nielsen, Torsten O.
#'   and Ellis, Matthew J. and Perou, Charles M. and Bernard, Philip S. (2009)
#'   "Supervised Risk Predictor of Breast Cancer Based on Intrinsic Subtypes",
#'   Journal of Clinical Oncology, 27(8):1160-1167
#' Desmedt C, Haibe-Kains B, Wirapati P, Buyse M, Larsimont D, Bontempi G, Delorenzi
#'   M, Piccart M, and Sotiriou C (2008) "Biological processes associated with breast
#'   cancer clinical outcome depend on the molecular subtypes", Clinical Cancer
#'   Research, 14(16):5158-5165.
#' Wirapati P, Sotiriou C, Kunkel S, Farmer P, Pradervand S, Haibe-Kains B, Desmedt
#'   C, Ignatiadis M, Sengstag T, Schutz F, Goldstein DR, Piccart MJ and Delorenzi M
#'   (2008) "Meta-analysis of Gene-Expression Profiles in Breast Cancer: Toward a
#'   Unified Understanding of Breast Cancer Sub-typing and Prognosis Signatures",
#'   Breast Cancer Research, 10(4):R65.
#' Haibe-Kains B, Desmedt C, Loi S, Culhane AC, Bontempi G, Quackenbush J, Sotiriou
#'   C. (2012) "A three-gene model to robustly identify breast cancer molecular
#'   subtypes.", J Natl Cancer Inst., 104(4):311-325.
#' Curtis C, Shah SP, Chin SF, Turashvili G, Rueda OM, Dunning MJ, Speed D, Lynch AG,
#'   Samarajiwa S, Yuan Y, Graf S, Ha G, Haffari G, Bashashati A, Russell R, McKinney
#'   S; METABRIC Group, Langerod A, Green A, Provenzano E, Wishart G, Pinder S, Watson
#'   P, Markowetz F, Murphy L, Ellis I, Purushotham A, Borresen-Dale AL, Brenton JD,
#'   Tavare S, Caldas C, Aparicio S. (2012) "The genomic and transcriptomic
#'   architecture of 2,000 breast tumours reveals novel subgroups.", Nature,
#'   486(7403):346-352.
#' Paquet ER, Hallett MT. (2015) "Absolute assignment of breast cancer intrinsic
#'   molecular subtype.", J Natl Cancer Inst., 107(1):357.
#' Aleix Prat, Joel S Parker, Olga Karginova, Cheng Fan, Chad Livasy, Jason I
#'   Herschkowitz, Xiaping He, and Charles M. Perou (2010) "Phenotypic and molecular
#'   characterization of the claudin-low intrinsic subtype of breast cancer", Breast
#'   Cancer Research, 12(5):R68
#'
#' @seealso
#' [genefu::subtype.cluster.predict], [genefu::intrinsic.cluster.predict]
#'
#' @examples
#' ##### without mapping (affy hgu133a or plus2 only)
#' # load VDX data
#' data(vdxs)
#' data(AIMSmodel)
#'
#' # Subtype Clustering Model fitted on EXPO and applied on VDX
#' sbt.vdx.SCMGENE <- molecular.subtyping(sbt.model="scmgene",
#'   data=data.vdxs, annot=annot.vdxs, do.mapping=FALSE)
#' table(sbt.vdx.SCMGENE$subtype)
#'
#' # Using the AIMS molecular subtyping algorithm
#' sbt.vdxs.AIMS <- molecular.subtyping(sbt.model="AIMS", data=data.vdxs,
#'                                      annot=annot.vdxs, do.mapping=FALSE)
#' table(sbt.vdxs.AIMS$subtype)
#'
#' # Using the IntClust molecular subtyping algorithm
#' colnames(annot.vdxs)[3]<-"Gene.Symbol"
#' sbt.vdxs.intClust <- molecular.subtyping(sbt.model="intClust", data=data.vdxs,
#'   annot=annot.vdxs, do.mapping=FALSE)
#' table(sbt.vdxs.intClust$subtype)
#'
#' ##### with mapping
#' # load NKI data
#' data(nkis)
#'
#' # Subtype Clustering Model fitted on EXPO and applied on NKI
#' sbt.nkis <- molecular.subtyping(sbt.model="scmgene", data=data.nkis,
#'   annot=annot.nkis, do.mapping=TRUE)
#' table(sbt.nkis$subtype)
#'
#' ##### with mapping
#' ## load vdxs data
#' data(vdxs)
#'
#' ## Claudin-Low classification of 150 VDXS samples
#' sbt.vdxs.CL <- molecular.subtyping(sbt.model="claudinLow", data=data.vdxs,
#'   annot=annot.vdxs, do.mapping=TRUE)
#' table(sbt.vdxs.CL$subtype)
#'
#' @md
#' @export
molecular.subtyping <- function(sbt.model=c("scmgene", "scmod1", "scmod2",
  "pam50", "ssp2006", "ssp2003", "intClust", "AIMS","claudinLow"), data, annot,
  do.mapping=FALSE, verbose=FALSE)
{

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
        sbts <- subtype.cluster.predict(sbt.model=scmgene.robust, data=data,
          annot=annot, do.mapping=do.mapping)[c("subtype2", "subtype.proba2")]
      },
      "scmod1" = {
        sbts <- subtype.cluster.predict(sbt.model=scmod1.robust, data=data,
          annot=annot, do.mapping=do.mapping)[c("subtype2", "subtype.proba2")]
      },
      "scmod2" = {
        sbts <- subtype.cluster.predict(sbt.model=scmod2.robust, data=data,
          annot=annot, do.mapping=do.mapping)[c("subtype2", "subtype.proba2")]
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
      names(gid1) <- paste("geneid", rownames(train$xd), sep=".")
      gid2 <- as.numeric(as.character(annot[ ,"EntrezGene.ID"]))
      names(gid2) <- colnames(data)

      ## remove missing and duplicated geneids from the gene list
      rm.ix <- is.na(gid1) | duplicated(gid1)
      gid1 <- gid1[!rm.ix]
      rr <- geneid.map(geneid1=gid2, data1=data, geneid2=gid1, verbose=FALSE)
      gt <- length(rr$geneid2)
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

    test <- medianCtr(t(data)) #probes as rows, median-centered
    #Run Classifier Call
	train2 <- train$xd
	rownames(train2) <- paste("geneid", rownames(train2), sep=".")
    predout <- claudinLow(x=train2, classes=as.matrix(train$classes$Group,ncol=1), y=test)
    sbts <- NULL
    sbts$subtype <- factor(as.character(predout$predictions$Call))
    colnames(predout$centroids) <- c("Claudin","Others")

    ## apply the nearest centroid classifier to classify the samples again
    ncor <- t(apply(X=data, MARGIN=1, FUN=function(x, y) {
      rr <- array(NA, dim=ncol(y), dimnames=list(colnames(y)))
      if (sum(complete.cases(x, y)) > 3) {
        rr <- cor(x=x, y=y, method="spearman", use="complete.obs")
      }
      return (rr)
    }, y=predout$centroids))

    #Calculate posterior probability based on the correlationss
   # nproba <- t(apply(X=ncor, MARGIN=1, FUN=function(x) { return(abs(x) / sum(abs(x), na.rm=TRUE)) }))

    # negative correlations are truncated to zero since they have no meaning for subtypes identification
    nproba <- t(apply(X=ncor, MARGIN=1, FUN=function (x) {
      rr <- array(NA, dim=length(x), dimnames=list(names(x)))
      x[!is.na(x) & x < 0] <- 0
      if (!all(is.na(x))) {
        rr <- x / sum(x, na.rm=TRUE)
      }
      return (rr)
    }))

    sbts$subtype.proba<-nproba

    ## compute crisp classification - in this case, really based on the binary call from the CL classifier
    #     sbts$subtype.crisp <- t(
    #       apply(sbts$subtype.proba, 1, function (x) {
    #         xx <- array(0, dim=length(x), dimnames=list(names(x)))
    #         xx[which.max(x)] <- 1
    #         return (xx)
    #       })
    #     )
    #     colnames(sbts$subtype.crisp)<-c("Claudin","Others")

    # In this case, really based on the binary call from the CL classifier. Use that for accuracy
      CLsubtypes<-c("Claudin","Others")
      sbts$subtype.crisp <- matrix(0, nrow=nrow(predout$predictions), ncol=2,dimnames=list(rownames(predout$predictions),CLsubtypes))
      for(count in 1:nrow(predout$predictions))
      {
        if(predout$predictions$Call[count]=="Others")
          sbts$subtype.crisp[count,2]<-1
        else sbts$subtype.crisp[count,1]<-1
      }
  }
  return (sbts)
}
