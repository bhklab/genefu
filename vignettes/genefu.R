## ----install, eval=FALSE, results='hide', message=FALSE-----------------------
#  BiocManager::install("genefu")

## ----load, eval=TRUE, results='hide', message=FALSE---------------------------
library(genefu)
library(xtable)
library(rmeta)
library(Biobase)
library(caret)

## ----install_data, eval=FALSE, results='hide', message=FALSE------------------
#  BiocManager::install("breastCancerMAINZ")
#  BiocManager::install("breastCancerTRANSBIG")
#  BiocManager::install("breastCancerUPP")
#  BiocManager::install("breastCancerUNT")
#  BiocManager::install("breastCancerNKI")

## ----load_data, eval=TRUE, results='hide', message=FALSE----------------------
library(breastCancerMAINZ)
library(breastCancerTRANSBIG)
library(breastCancerUPP)
library(breastCancerUNT)
library(breastCancerNKI)

## ----findDuplicatedPatients,eval=TRUE,results='hide',message=FALSE------------
data(breastCancerData)
cinfo <- colnames(pData(mainz7g))
data.all <- c("transbig7g"=transbig7g, "unt7g"=unt7g, "upp7g"=upp7g,
              "mainz7g"=mainz7g, "nki7g"=nki7g)

idtoremove.all <- NULL
duplres <- NULL

## No overlaps in the MainZ and NKI datasets.

## Focus on UNT vs UPP vs TRANSBIG
demo.all <- rbind(pData(transbig7g), pData(unt7g), pData(upp7g))
dn2 <- c("TRANSBIG", "UNT", "UPP")

## Karolinska
## Search for the VDXKIU, KIU, UPPU series
ds2 <- c("VDXKIU", "KIU", "UPPU")
demot <- demo.all[complete.cases(demo.all[ , c("series")]) &
                    is.element(demo.all[ , "series"], ds2), ]

# Find the duplicated patients in that series
duplid <- sort(unique(demot[duplicated(demot[ , "id"]), "id"]))
duplrest <- NULL
for(i in 1:length(duplid)) {
  tt <- NULL
  for(k in 1:length(dn2)) {
    myx <- sort(row.names(demot)[complete.cases(demot[ , c("id", "dataset")]) &
                                   demot[ , "id"] == duplid[i] & 
                                   demot[ , "dataset"] == dn2[k]])
    if(length(myx) > 0) { tt <- c(tt, myx) }
  }
  duplrest <- c(duplrest, list(tt))
}
names(duplrest) <- duplid
duplres <- c(duplres, duplrest)

## Oxford
## Search for the VVDXOXFU, OXFU series
ds2 <- c("VDXOXFU", "OXFU")
demot <- demo.all[complete.cases(demo.all[ , c("series")]) & 
                    is.element(demo.all[ , "series"], ds2), ]

# Find the duplicated patients in that series
duplid <- sort(unique(demot[duplicated(demot[ , "id"]), "id"]))
duplrest <- NULL
for(i in 1:length(duplid)) {
  tt <- NULL
  for(k in 1:length(dn2)) {
    myx <- sort(row.names(demot)[complete.cases(demot[ , c("id", "dataset")]) &
                                   demot[ , "id"] == duplid[i] & 
                                   demot[ , "dataset"] == dn2[k]])
    if(length(myx) > 0) { tt <- c(tt, myx) }
  }
  duplrest <- c(duplrest, list(tt))
}
names(duplrest) <- duplid
duplres <- c(duplres, duplrest)

## Full set duplicated patients
duPL <- sort(unlist(lapply(duplres, function(x) { return(x[-1]) } )))

## ----CalculateMolecularSubtypes-----------------------------------------------
dn <- c("transbig", "unt", "upp", "mainz", "nki")
dn.platform <- c("affy", "affy", "affy", "affy", "agilent")
res <- ddemo.all <- ddemo.coln <- NULL

for(i in 1:length(dn)) {

  ## load dataset
  dd <- get(data(list=dn[i]))
  #Remove duplicates identified first
  message("obtained dataset!")

  #Extract expression set, pData, fData for each dataset
  ddata <- t(exprs(dd))

  ddemo <- phenoData(dd)@data

  if(length(intersect(rownames(ddata),duPL))>0)
  {
  ddata<-ddata[-which(rownames(ddata) %in% duPL),]
  ddemo<-ddemo[-which(rownames(ddemo) %in% duPL),]
  }

  dannot <- featureData(dd)@data

  # MOLECULAR SUBTYPING
  # Perform subtyping using scmod2.robust
  # scmod2.robust: List of parameters defining the subtype clustering model
  # (as defined by Wirapati et al)

  # OBSOLETE FUNCTION CALL - OLDER VERSIONS OF GENEFU
  # SubtypePredictions<-subtype.cluster.predict(sbt.model=scmod2.robust,data=ddata,
  #                                               annot=dannot,do.mapping=TRUE,
  #                                               verbose=TRUE)

  # CURRENT FUNCTION CALL - NEWEST VERSION OF GENEFU
  SubtypePredictions <- molecular.subtyping(sbt.model = "scmod2",data = ddata,
                                          annot = dannot,do.mapping = TRUE)

  #Get sample counts pertaining to each subtype
  table(SubtypePredictions$subtype)
  #Select samples pertaining to Basal Subtype
  Basals<-names(which(SubtypePredictions$subtype == "ER-/HER2-"))
  #Select samples pertaining to HER2 Subtype
  HER2s<-names(which(SubtypePredictions$subtype == "HER2+"))
  #Select samples pertaining to Luminal Subtypes
  LuminalB<-names(which(SubtypePredictions$subtype == "ER+/HER2- High Prolif"))
  LuminalA<-names(which(SubtypePredictions$subtype == "ER+/HER2- Low Prolif"))

  #ASSIGN SUBTYPES TO EVERY SAMPLE, ADD TO THE EXISTING PHENODATA
  ddemo$SCMOD2<-SubtypePredictions$subtype
  ddemo[LuminalB,]$SCMOD2<-"LumB"
  ddemo[LuminalA,]$SCMOD2<-"LumA"
  ddemo[Basals,]$SCMOD2<-"Basal"
  ddemo[HER2s,]$SCMOD2<-"Her2"

  # Perform subtyping using PAM50
  # Matrix should have samples as ROWS, genes as COLUMNS
  # rownames(dannot)<-dannot$probe<-dannot$EntrezGene.ID

  # OLDER FUNCTION CALL
  # PAM50Preds<-intrinsic.cluster.predict(sbt.model=pam50,data=ddata,
  #                                         annot=dannot,do.mapping=TRUE,
  #                                         verbose=TRUE)

  # NEWER FUNCTION CALL BASED ON MOST RECENT VERSION
  PAM50Preds<-molecular.subtyping(sbt.model = "pam50",data=ddata,
                                        annot=dannot,do.mapping=TRUE)


  table(PAM50Preds$subtype)
  ddemo$PAM50<-PAM50Preds$subtype
  LumA<-names(PAM50Preds$subtype)[which(PAM50Preds$subtype == "LumA")]
  LumB<-names(PAM50Preds$subtype)[which(PAM50Preds$subtype == "LumB")]
  ddemo[LumA,]$PAM50<-"LumA"
  ddemo[LumB,]$PAM50<-"LumB"

  ddemo.all <- rbind(ddemo, ddemo.all)
}

## ----CompareMolecularSubtypesByConfusionMatrix--------------------------------
# Obtain the subtype prediction counts for PAM50
table(ddemo.all$PAM50)
Normals<-rownames(ddemo.all[which(ddemo.all$PAM50 == "Normal"),])

# Obtain the subtype prediction counts for SCMOD2
table(ddemo.all$SCMOD2)

ddemo.all$PAM50<-as.character(ddemo.all$PAM50)
# We compare the samples that are predicted as pertaining to a molecular subtyp
# We ignore for now the samples that predict as 'Normal' by PAM50
confusionMatrix(
  factor(ddemo.all[-which(rownames(ddemo.all) %in% Normals),]$SCMOD2),
  factor(ddemo.all[-which(rownames(ddemo.all) %in% Normals),]$PAM50)
  )

## ----CompareSurvivalBySubtypes------------------------------------------------
# http://www.inside-r.org/r-doc/survival/survfit.coxph
library(survival)
ddemo<-ddemo.all
data.for.survival.SCMOD2 <- ddemo[,c("e.os", "t.os", "SCMOD2","age")]
data.for.survival.PAM50 <- ddemo[,c("e.os", "t.os", "PAM50","age")]
# Remove patients with missing survival information
data.for.survival.SCMOD2 <- 
  data.for.survival.SCMOD2[complete.cases(data.for.survival.SCMOD2),]
data.for.survival.PAM50 <- 
  data.for.survival.PAM50[complete.cases(data.for.survival.PAM50),]

days.per.month <- 30.4368
days.per.year <- 365.242

data.for.survival.PAM50$months_to_death <- 
  data.for.survival.PAM50$t.os / days.per.month
data.for.survival.PAM50$vital_status <- data.for.survival.PAM50$e.os == "1"
surv.obj.PAM50 <- survfit(Surv(data.for.survival.PAM50$months_to_death,
                               data.for.survival.PAM50$vital_status) ~ 
                            data.for.survival.PAM50$PAM50)

data.for.survival.SCMOD2$months_to_death <- 
  data.for.survival.SCMOD2$t.os / days.per.month
data.for.survival.SCMOD2$vital_status <- data.for.survival.SCMOD2$e.os == "1"
surv.obj.SCMOD2 <- survfit(Surv(
  data.for.survival.SCMOD2$months_to_death,
  data.for.survival.SCMOD2$vital_status) ~ data.for.survival.SCMOD2$SCMOD2)

message("KAPLAN-MEIR CURVE - USING PAM50")

plot(main = "Surival Curves PAM50", surv.obj.PAM50,
     col =c("#006d2c", "#8856a7","#a50f15", "#08519c", "#000000"),lty = 1,lwd = 3,
     xlab = "Time (months)",ylab = "Probability of Survival")
legend("topright",
       fill = c("#006d2c", "#8856a7","#a50f15", "#08519c", "#000000"),
       legend = c("Basal","Her2","LumA","LumB","Normal"),bty = "n")

message("KAPLAN-MEIR CURVE - USING SCMOD2")

plot(main = "Surival Curves SCMOD2", surv.obj.SCMOD2,
     col =c("#006d2c", "#8856a7","#a50f15", "#08519c"),lty = 1,lwd = 3,
     xlab = "Time (months)",ylab = "Probability of Survival")
legend("topright",
       fill = c("#006d2c", "#8856a7","#a50f15", "#08519c"),
       legend = c("Basal","Her2","LumA","LumB"),bty = "n")

## GENERATE A OVERLAYED PLOT OF SURVIVAL CURVES
message("Overlayed Surival Plots based on PAM50 and SCMOD2")
                          ## Basal    Her2        LuminalA  LuminalB   Normal
plot(surv.obj.PAM50,col =c("#006d2c", "#8856a7","#a50f15", "#08519c", "#000000"),lty = 1,lwd = 3,
     xlab = "Time (months)",ylab = "Probability of Survival",ymin = 0.2)
legend("topright",
       fill = c("#006d2c", "#8856a7","#a50f15", "#08519c", "#000000"),
       legend = c("Basal","Her2","LumA","LumB","Normal"),bty = "n")

par(new=TRUE)
                            ## Basal    Her2        LuminalA  LuminalB
lines(surv.obj.SCMOD2,col =c("#006d2c", "#8856a7","#a50f15", "#08519c"),lwd=2,lty=5)
legend("bottomright",c("PAM50","SCMOD2"),lty=c("solid", "dashed"))

## ----CalculatedCVPL-----------------------------------------------------------
set.seed(12345)

PAM5_CVPL<-cvpl(x=data.for.survival.PAM50$age,
                surv.time=data.for.survival.PAM50$months_to_death,
                surv.event=data.for.survival.PAM50$vital_status,
                strata=as.integer(factor(data.for.survival.PAM50$PAM50)),
                nfold=10, setseed=54321)$cvpl

SCMOD2_CVPL<-cvpl(x=data.for.survival.SCMOD2$age,
                    surv.time=data.for.survival.SCMOD2$months_to_death,
                    surv.event=data.for.survival.SCMOD2$vital_status,
                    strata=as.integer(factor(data.for.survival.SCMOD2$SCMOD2)),
                    nfold=10, setseed=54321)$cvpl

print.data.frame(data.frame(cbind(PAM5_CVPL,SCMOD2_CVPL)))

## ----computeRiskScore---------------------------------------------------------
dn <- c("transbig", "unt", "upp", "mainz", "nki")
dn.platform <- c("affy", "affy", "affy", "affy", "agilent")

res <- ddemo.all <- ddemo.coln <- NULL
for(i in 1:length(dn)) {

  ## load dataset
  dd <- get(data(list=dn[i]))

  #Extract expression set, pData, fData for each dataset
  ddata <- t(exprs(dd))
  ddemo <- phenoData(dd)@data
  dannot <- featureData(dd)@data
  ddemo.all <- c(ddemo.all, list(ddemo))
  if(is.null(ddemo.coln))
  { ddemo.coln <- colnames(ddemo) } else
  { ddemo.coln <- intersect(ddemo.coln, colnames(ddemo)) }
  rest <- NULL

  ## AURKA
  ## if affy platform consider the probe published in Desmedt et al., CCR, 2008
  if(dn.platform[i] == "affy") { domap <- FALSE } else { domap <- TRUE }
  modt <- scmgene.robust$mod$AURKA
  ## if agilent platform consider the probe published in Desmedt et al., CCR, 2008
  if(dn.platform[i] == "agilent") {
    domap <- FALSE
    modt[ , "probe"] <- "NM_003600"
  }
  rest <- cbind(rest, "AURKA"=sig.score(x=modt, data=ddata, annot=dannot, 
                                        do.mapping=domap)$score)

  ## ESR1
  ## if affy platform consider the probe published in Desmedt et al., CCR, 2008
  if(dn.platform[i] == "affy") { domap <- FALSE } else { domap <- TRUE }
  modt <- scmgene.robust$mod$ESR1
  ## if agilent platform consider the probe published in Desmedt et al., CCR, 2008
  if(dn.platform[i] == "agilent") {
    domap <- FALSE
    modt[ , "probe"] <- "NM_000125"
  }
  rest <- cbind(rest, "ESR1"=sig.score(x=modt, data=ddata, annot=dannot, 
                                       do.mapping=domap)$score)

  ## ERBB2
  ## if affy platform consider the probe published in Desmedt et al., CCR, 2008
  if(dn.platform[i] == "affy") { domap <- FALSE } else { domap <- TRUE }
  modt <- scmgene.robust$mod$ERBB2
  ## if agilent platform consider the probe published in Desmedt et al., CCR, 2008
  if(dn.platform[i] == "agilent") {
    domap <- FALSE
    modt[ , "probe"] <- "NM_004448"
  }
  rest <- cbind(rest, "ERBB2"=sig.score(x=modt, data=ddata, annot=dannot, 
                                        do.mapping=domap)$score)

  ## NPI
  ss <- ddemo[ , "size"]
  gg <- ddemo[ , "grade"]
  nn <- rep(NA, nrow(ddemo))
  nn[complete.cases(ddemo[ , "node"]) & ddemo[ , "node"] == 0] <- 1
  nn[complete.cases(ddemo[ , "node"]) & ddemo[ , "node"] == 1] <- 3
  names(ss) <- names(gg) <- names(nn) <- rownames(ddemo)
  rest <- cbind(rest, "NPI"=npi(size=ss, grade=gg, node=nn, na.rm=TRUE)$score)

  ## GGI
  if(dn.platform[i] == "affy") { domap <- FALSE } else { domap <- TRUE }
  rest <- cbind(rest, "GGI"=ggi(data=ddata, annot=dannot, 
                                do.mapping=domap)$score)

  ## GENIUS
  if(dn.platform[i] == "affy") { domap <- FALSE } else { domap <- TRUE }
  rest <- cbind(rest, "GENIUS"=genius(data=ddata, annot=dannot, 
                                      do.mapping=domap)$score)

  ## ENDOPREDICT
  if(dn.platform[i] == "affy") { domap <- FALSE } else { domap <- TRUE }
  rest <- cbind(rest, "EndoPredict"=endoPredict(data=ddata, annot=dannot, 
                                                do.mapping=domap)$score)

  # OncotypeDx
  if(dn.platform[i] == "affy") { domap <- FALSE } else { domap <- TRUE }
  rest <- cbind(rest, "OncotypeDx"=oncotypedx(data=ddata, annot=dannot, 
                                              do.mapping=domap)$score)

  ## TamR
  # Note: risk is not implemented, the function will return NA values
  if(dn.platform[i] == "affy") { domap <- FALSE } else { domap <- TRUE }
  rest <- cbind(rest, "TAMR13"=tamr13(data=ddata, annot=dannot, 
                                      do.mapping=domap)$score)

  ## GENE70
  # Need to do mapping for Affy platforms because this is based on Agilent.
  # Hence the mapping rule is reversed here!
  if(dn.platform[i] == "affy") { domap <- TRUE } else { domap <- FALSE }
  rest <- cbind(rest, "GENE70"=gene70(data=ddata, annot=dannot, std="none",
                                      do.mapping=domap)$score)

  ## Pik3cags
  if(dn.platform[i] == "affy") { domap <- FALSE } else { domap <- TRUE }
  rest <- cbind(rest, "PIK3CA"=pik3cags(data=ddata, annot=dannot, 
                                        do.mapping=domap))

  ## rorS
  # Uses the pam50 algorithm. Need to do mapping for both Affy and Agilent
  rest <- cbind(rest, "rorS"=rorS(data=ddata, annot=dannot, 
                                  do.mapping=TRUE)$score)

  ## GENE76
  # Mainly designed for Affy platforms. Has been excluded here

  # BIND ALL TOGETHER
  res <- rbind(res, rest)
}
names(ddemo.all) <- dn

## ----simplifyAndRemoveDuplicatePatients---------------------------------------
ddemot <- NULL
for(i in 1:length(ddemo.all)) {
  ddemot <- rbind(ddemot, ddemo.all[[i]][ , ddemo.coln, drop=FALSE])
}
res[complete.cases(ddemot[ ,"dataset"]) & ddemot[ ,"dataset"] == "VDX", "GENIUS"] <- NA

## select only untreated node-negative patients with all risk predictions
## ie(incomplete cases (where risk prediction may be missing for a sample) are subsequently removed))
# Note that increasing the number of risk prediction analyses
# may increase the number of incomplete cases
# In the previous vignette for genefu version1, we were only testing 4 risk predictors,
# so we had a total of 722 complete cases remaining
# Here, we are now testing 12 risk predictors, so we only have 713 complete cases remaining.
# The difference of 9 cases between the two versions are all from the NKI dataset.
myx <- complete.cases(res, ddemot[ , c("node", "treatment")]) &
  ddemot[ , "treatment"] == 0 & ddemot[ , "node"] == 0 & !is.element(rownames(ddemot), duPL)

res <- res[myx, , drop=FALSE]
ddemot <- ddemot[myx, , drop=FALSE]

## ----cindexComputation--------------------------------------------------------
cc.res <- complete.cases(res)
datasetList <- c("MAINZ","TRANSBIG","UPP","UNT","NKI")
riskPList <- c("AURKA","ESR1","ERBB2","NPI", "GGI", "GENIUS",
               "EndoPredict","OncotypeDx","TAMR13","GENE70","PIK3CA","rorS")
setT <- setE <- NULL
resMatrix <- as.list(NULL)

for(i in datasetList)
{
  dataset.only <- ddemot[,"dataset"] == i
  patientsAll <- cc.res & dataset.only

  ## set type of available survival data
  if(i == "UPP") {
    setT <- "t.rfs"
    setE <- "e.rfs"
  } else {
    setT <- "t.dmfs"
    setE <- "e.dmfs"
  }

  # Calculate cindex computation for each predictor
  for (Dat in riskPList)
  {
    cindex <- t(apply(X=t(res[patientsAll,Dat]), MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
    return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=ddemot[patientsAll,setT], z=ddemot[patientsAll, setE]))

    resMatrix[[Dat]] <- rbind(resMatrix[[Dat]], cindex)
  }
}

## ----combineEstimations-------------------------------------------------------
for(i in names(resMatrix)){
  #Get a meta-estimate
  ceData <- combine.est(x=resMatrix[[i]][,"cindex"], x.se=resMatrix[[i]][,"cindex.se"], hetero=TRUE)
  cLower <- ceData$estimate + qnorm(0.025, lower.tail=TRUE) * ceData$se
  cUpper <- ceData$estimate + qnorm(0.025, lower.tail=FALSE) * ceData$se

  cindexO <- cbind("cindex"=ceData$estimate, "cindex.se"=ceData$se, "lower"=cLower, "upper"=cUpper)
  resMatrix[[i]] <- rbind(resMatrix[[i]], cindexO)
  rownames(resMatrix[[i]]) <- c(datasetList, "Overall")
}

## ----computePValues-----------------------------------------------------------
pv <- sapply(resMatrix, function(x) { return(x["Overall", c("cindex","cindex.se")]) })
pv <- apply(pv, 2, function(x) { return(pnorm((x[1] - 0.5) / x[2], lower.tail=x[1] < 0.5)) })
printPV <- matrix(pv,ncol=length(names(resMatrix)))
rownames(printPV) <- "P-value"
colnames(printPV) <- names(pv)
printPV<-t(printPV)

## ----printPvalue,results="asis"-----------------------------------------------
xtable(printPV, digits=c(0, -1))

## ----forestplotDatasets,echo=TRUE---------------------------------------------
RiskPList <- c("AURKA","ESR1","ERBB2","NPI", "GGI", "GENIUS",
               "EndoPredict","OncotypeDx","TAMR13","GENE70","PIK3CA","rorS")
datasetListF <- c("MAINZ","TRANSBIG","UPP","UNT","NKI", "Overall")
myspace <- "   "
par(mfrow=c(2,2))
  for (RP in RiskPList)
  {

  #<<forestplotDat,fig=TRUE>>=
  ## Forestplot
  tt <- rbind(resMatrix[[RP]][1:5,],
            "Overall"=resMatrix[[RP]][6,])

  tt <- as.data.frame(tt)
  labeltext <- (datasetListF)

  r.mean <- c(tt$cindex)
  r.lower <- c(tt$lower)
  r.upper <- c(tt$upper)

  metaplot.surv(mn=r.mean, lower=r.lower, upper=r.upper, labels=labeltext, xlim=c(0.3,0.9),
                boxsize=0.5, zero=0.5,
                col=meta.colors(box="royalblue",line="darkblue",zero="firebrick"),
                main=paste(RP))

  }

## ----forestplotOverall,echo=TRUE----------------------------------------------
## Overall Forestplot
mybigspace <- "       "
tt <- rbind("OverallA"=resMatrix[["AURKA"]][6,],
            "OverallE1"=resMatrix[["ESR1"]][6,],
            "OverallE2"=resMatrix[["ERBB2"]][6,],
            "OverallN"=resMatrix[["NPI"]][6,],
          "OverallM"=resMatrix[["GGI"]][6,],
          "OverallG"=resMatrix[["GENIUS"]][6,],
          "OverallE3"=resMatrix[["EndoPredict"]][6,],
          "OverallOD"=resMatrix[["OncotypeDx"]][6,],
          "OverallT"=resMatrix[["TAMR13"]][6,],
          "OverallG70"=resMatrix[["GENE70"]][6,],
          "OverallP"=resMatrix[["PIK3CA"]][6,],
          "OverallR"=resMatrix[["rorS"]][6,]
          )

tt <- as.data.frame(tt)
labeltext <- cbind(c("Risk Prediction","AURKA","ESR1","ERBB2","NPI",
                     "GGI","GENIUS","EndoPredict","OncotypeDx","TAMR13","GENE70","PIK3CA","rorS"))

r.mean <- c(NA,tt$cindex)
r.lower <- c(NA,tt$lower)
r.upper <- c(NA,tt$upper)

metaplot.surv(mn=r.mean, lower=r.lower, upper=r.upper, labels=labeltext, xlim=c(0.35,0.75),
              boxsize=0.5, zero=0.5,
              col=meta.colors(box="royalblue",line="darkblue",zero="firebrick"),
              main="Overall Concordance Index")

## ----computeCindexWithPvalue--------------------------------------------------
cc.res <- complete.cases(res)
datasetList <- c("MAINZ","TRANSBIG","UPP","UNT","NKI")
riskPList <- c("AURKA","ESR1","ERBB2","NPI","GGI","GENIUS",
               "EndoPredict","OncotypeDx","TAMR13","GENE70","PIK3CA","rorS")
setT <- setE <- NULL
resMatrixFull <- as.list(NULL)

for(i in datasetList)
{
  dataset.only <- ddemot[,"dataset"] == i
  patientsAll <- cc.res & dataset.only

  ## set type of available survival data
  if(i == "UPP") {
    setT <- "t.rfs"
    setE <- "e.rfs"
  } else {
    setT <- "t.dmfs"
    setE <- "e.dmfs"
  }

  ## cindex and p-value computation per algorithm
  for (Dat in riskPList)
  {
    cindex <- t(apply(X=t(res[patientsAll,Dat]), MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
    return(tt); },
    y=ddemot[patientsAll,setT], z=ddemot[patientsAll, setE]))

    resMatrixFull[[Dat]] <- rbind(resMatrixFull[[Dat]], cindex)
  }
}

for(i in names(resMatrixFull)){
  rownames(resMatrixFull[[i]]) <- datasetList
}

ccmData <- tt <- rr <- NULL
for(i in 1:length(resMatrixFull)){
  tt <- NULL
  for(j in 1:length(resMatrixFull)){
    if(i != j) { rr <- cindex.comp.meta(list.cindex1=resMatrixFull[[i]],
                                        list.cindex2=resMatrixFull[[j]], hetero=TRUE)$p.value }
    else { rr <- 1 }
    tt <- cbind(tt, rr)
  }
  ccmData <- rbind(ccmData, tt)
}
ccmData <- as.data.frame(ccmData)
colnames(ccmData) <- riskPList
rownames(ccmData) <- riskPList

## ----computeCCMPval-----------------------------------------------------------
ccmDataPval <- matrix(p.adjust(data.matrix(ccmData), method="holm"),
                      ncol=length(riskPList), dimnames=list(rownames(ccmData),
                                                            colnames(ccmData)))

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

