#' @title Claudin-low classification for Breast Cancer Data
#'
#' @description
#' Subtyping method for identifying Claudin-Low Breast Cancer Samples.
#'   Code generously provided by Aleix Prat.
#'
#' @usage
#' claudinLow(x, classes="", y, nGenes="", priors="equal",
#'   std=FALSE, distm="euclidean", centroids=FALSE)
#'
#' @param x the data matrix of training samples, or pre-calculated centroids.
#' @param classes a list labels for use in coloring the points.
#' @param y the data matrix of test samples.
#' @param nGenes the number of genes selected when training the model.
#' @param priors 'equal' assumes equal class priors, 'class' calculates them
#'   based on proportion in the data.
#' @param std when true, the training and testing samples are standardized
#'   to mean=0 and var=1.
#' @param distm the distance metric for determining the nearest centroid,
#'   can be one of euclidean, pearson, or spearman.
#' @param centroids when true, it is assumed that x consists of pre-calculated centroids.
#'
#' @return
#' A list with items:
#' - predictions
#' - testData
#' - distances
#' - centroids
#'
#' @references
#' Aleix Prat, Joel S Parker, Olga Karginova, Cheng Fan, Chad Livasy, Jason
#'   I Herschkowitz, Xiaping He, and Charles M. Perou (2010) "Phenotypic and
#'   molecular characterization of the claudin-low intrinsic subtype of
#'   breast cancer", Breast Cancer Research, 12(5):R68
#'
#' @seealso
#' [genefu::medianCtr()], [genefu::claudinLowData]
#'
#' @examples
#' data(claudinLowData)
#'
#' #Training Set
#' train <- claudinLowData
#' train$xd <-  medianCtr(train$xd)
#' # Testing Set
#' test <- claudinLowData
#' test$xd <-  medianCtr(test$xd)
#'
#' # Generate Predictions
#' predout <- claudinLow(x=train$xd, classes=as.matrix(train$classes$Group,ncol=1), y=test$xd)
#'
#' # Obtain results
#' results <- cbind(predout$predictions, predout$distances)
#' #write.table(results,"T.E.9CELL.LINE_results.txt",sep="\t",col=T, row=FALSE)
#'
#' @md
#' @import limma stats utils
#' @export
claudinLow <- function(x, classes="", y, nGenes="", priors="equal", std=FALSE,
    distm="euclidean", centroids=FALSE){

  dataMatrix <- x
  features <- dim(x)[1]
  samples <- dim(x)[2]
  sampleNames <- dimnames(x)[[2]]
  featureNames <- dimnames(x)[[1]]

  #parse the test file - same as train file but no rows of classes
  tdataMatrix <- y
  tfeatures <- dim(y)[1]
  tsamples <- dim(y)[2]
  tsampleNames <- dimnames(y)[[2]]
  tfeatureNames <- dimnames(y)[[1]]

  #dimnames(tdataMatrix)[[2]] <- paste("x",seq(1,471))
  temp <- overlapSets(dataMatrix,tdataMatrix)
  dataMatrix <- temp$x
  tdataMatrix <- temp$y
  sfeatureNames <- row.names(dataMatrix)

  # standardize both sets
  if(std){
    dataMatrix <- standardize(dataMatrix)
    tdataMatrix <- standardize(tdataMatrix)
  }

  if(!centroids){
    thisClass <- as.vector(classes[,1])
    nClasses <- nlevels(as.factor(thisClass))
    classLevels <- levels(as.factor(thisClass))
    for(j in 1:nClasses){
      thisClass[thisClass==classLevels[j]] <- j
    }
    thisClass <- as.numeric(thisClass)
    dataMatrix <- dataMatrix[,!(is.na(thisClass))]
    thisClass <- thisClass[!(is.na(thisClass))]

    scores <- apply(dataMatrix,1,limma::bwss,thisClass)
    trainscores <- vector()
    for(j in 1:dim(dataMatrix)[1]){
      trainscores[j] <- scores[[row.names(dataMatrix)[j]]]$bss / scores[[row.names(dataMatrix)[j]]]$wss
    }

    dataMatrix <- dataMatrix[sort.list(trainscores,decreasing=TRUE),]
    tdataMatrix <- tdataMatrix[sort.list(trainscores,decreasing=TRUE),]

    if(nGenes==""){
      nGenes <- dim(dataMatrix)[1]
    }
    print(paste("Number of genes used:",nGenes))

    dataMatrix <- dataMatrix[1:nGenes,]
    tdataMatrix <- tdataMatrix[1:nGenes,]

    centroids <- matrix(,nrow=nGenes,ncol=nClasses)
    for(j in 1:nClasses){
      centroids[,j] <- apply(dataMatrix[,thisClass==j],1,mean)
    }
    dimnames(centroids) <- list(row.names(dataMatrix),NULL)

  }else{
    nGenes <- dim(dataMatrix)[1]
    print(paste("Number of genes used:",nGenes))
    centroids <- dataMatrix
    nClasses <- dim(centroids)[2]
    classLevels <- dimnames(centroids)[[2]]
  }

  distances <- matrix(ncol=nClasses,nrow=dim(tdataMatrix)[2])
  for(j in 1:nClasses){
    if(distm=="euclidean"){
      distances[,j] <- dist(t(cbind(centroids[,j],tdataMatrix)))[1:(dim(tdataMatrix)[2])]
    }
    if(distm=="correlation" | distm=="pearson"){
      distances[,j] <- t(-1*cor(cbind(centroids[,j],tdataMatrix),use="pairwise.complete.obs"))[2:(dim(tdataMatrix)[2]+1)]
    }
    if(distm=="spearman"){
      distances[,j] <- t(-1*cor(cbind(centroids[,j],tdataMatrix),method="spearman",use="pairwise.complete.obs"))[2:(dim(tdataMatrix)[2]+1)]
    }
    colnames(distances) <- c("euclidian distance to Claudin-low", "euclidian distance to Others")
    rownames(distances) <- tsampleNames

  }

  scores <- apply(distances,1,min)
  prediction <- vector(length=tsamples)
  for(i in 1:tsamples){
    prediction[i] <- classLevels[match(scores[i],distances[i,])]
  }
  names(prediction) <- tsampleNames
  prediction <- data.frame(Samples=tsampleNames, prediction)
  colnames(prediction) <- c("Samples", "Call")
  return(list(predictions=prediction,testData=tdataMatrix,distances=distances,centroids=centroids))
}

