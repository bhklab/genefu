#' @title Function to identify bimodality for gene expression or signature score
#'
#' @description
#' This function fits a mixture of two Gaussians to identify bimodality. 
#' Useful to identify ER of HER2 status of breast tumors using 
#' ESR1 and ERBB2 expressions respectively.
#'
#' @usage
#' bimod(x, data, annot, do.mapping = FALSE, mapping, model = c("E", "V"), 
#' do.scale = TRUE, verbose = FALSE, ...)
#'   
#' @param x Matrix containing the gene(s) in the gene list in rows and at least three columns: 
#'   "probe", "EntrezGene.ID" and "coefficient" standing for the name of the probe, 
#'   the NCBI Entrez Gene id and the coefficient giving the direction and the strength 
#'   of the association of each gene in the gene list.
#' @param data Matrix of gene expressions with samples in rows and probes in columns, 
#'   dimnames being properly defined.
#' @param annot Matrix of annotations with at least one column named "EntrezGene.ID", 
#'   dimnames being properly defined.
#' @param do.mapping TRUE if the mapping through Entrez Gene ids must be performed (in case of 
#'   ambiguities, the most variant probe is kept for each gene), FALSE otherwise.
#' @param mapping Matrix with columns "EntrezGene.ID" and "probe" used to force the 
#'   mapping such that the probes are not selected based on their variance.
#' @param model Model name used in Mclust.
#' @param do.scale TRUE if the gene expressions or signature scores must be rescaled (see rescale), FALSE otherwise.
#' @param verbose TRUE to print informative messages, FALSE otherwise.
#' @param ... Additional parameters to pass to sig.score. 
#'
#' @return
#' A list with items:
#' - status: Status being 0 or 1.
#' - status1.proba: Probability p to be of status 1, the probability to 
#'   be of status 0 being 1-p.
#' - gaussians: Matrix of parameters fitted in the mixture of two 
#'   Gaussians. Matrix of NA values if EM algorithm did not converge.
#' - BIC: Values (gene expressions or signature scores) used to identify bimodality.
#' - BI: Bimodality Index (BI) as defined by Wang et al., 2009.
#' - x: Values (gene expressions or signature scores) used to identify bimodality
#'
#' @references
#' Desmedt C, Haibe-Kains B, Wirapati P, Buyse M, Larsimont D, Bontempi G, Delorenzi M, Piccart M,
#'   and Sotiriou C (2008) "Biological processes associated with breast cancer clinical outcome depend
#'   on the molecular subtypes", Clinical Cancer Research, 14(16):5158–5165.
#' Wirapati P, Sotiriou C, Kunkel S, Farmer P, Pradervand S, Haibe-Kains B, Desmedt C, Ignatiadis M,
#'   Sengstag T, Schutz F, Goldstein DR, Piccart MJ and Delorenzi M (2008) "Meta-analysis of 
#'   Gene-Expression Profiles in Breast Cancer: Toward a Unified Understanding of Breast Cancer Sub-typing
#'   and Prognosis Signatures", Breast Cancer Research, 10(4):R65.
#' Fraley C and Raftery E (2002) "Model-Based Clustering, Discriminant Analysis, and Density Estimation", 
#'   Journal of American Statistical Asscoiation, 97(458):611–631.
#' Wang J, Wen S, Symmans FW, Pusztai L and Coombes KR (2009) "The bimodality index: a criterion for 
#'   discovering and ranking bimodal signatures from cancer gene expression profiling data", Cancer 
#'   Informatics, 7:199–216.
#'   
#' @seealso
#' [mclust::Mclust]
#' 
#' @examples
#' # load NKI data
#' data(nkis)
#' # load gene modules from Desmedt et al. 2008
#' data(mod1)
#' # retrieve esr1 affy probe and Entrez Gene id 
#' esr1 <- mod1$ESR1[1, ,drop=FALSE]
#' # computation of signature scores
#' esr1.bimod <- bimod(x=esr1, data=data.nkis, annot=annot.nkis, do.mapping=TRUE, 
#'   model="V", verbose=TRUE)
#' table("ER.IHC"=demo.nkis[ ,"er"], "ER.GE"=esr1.bimod$status)
#'
#' @md
#' @export
bimod <-
  function(x, data, annot, do.mapping=FALSE, mapping, model=c("E", "V"), do.scale=TRUE, verbose=FALSE, ...) {
    #require(mclust)
    model <- match.arg(model)
    dd <- sig.score(x=x, data=data, annot=annot, do.mapping=do.mapping, mapping=mapping, verbose=verbose, ...)$score
    if(do.scale) { dd <- (rescale(x=dd, q=0.05, na.rm=TRUE) - 0.5) * 2 }
    cc.ix <- complete.cases(dd)
    
    mystatus <- mystatus.proba <- rep(NA, nrow(data))
    names(mystatus) <- names(mystatus.proba) <- dimnames(data)[[1]]
    res <- matrix(NA, nrow=3, ncol=2, dimnames=list(c("mean", "variance", "proportion"), paste("cluster", 1:2, sep=".")))
    mybic <- matrix(NA, nrow=10, ncol=1, dimnames=list(1:10, model))
    
    if(sum(cc.ix) >= 10) {	
      #How many Gaussians?
      rr <- mclust::Mclust(data=dd[cc.ix], modelNames=model, G=1:10)
      oo <- order(rr$BIC, decreasing=TRUE)[1]
      if(oo != 2) { warning(sprintf("%i is the most likely number of Gaussians!", oo)) }
      mybic <- rr$BIC
      
      #Only 2 Gaussians
      rr2 <- mclust::Mclust(data=dd[cc.ix], modelNames=model, G=2)
      if(is.null(rr2[[1]])) { ## EM algorithm did not converge
        return(list("status"=mystatus, "status1.proba"=mystatus.proba, "gaussians"=res, "BIC"=rr$BIC, "x"=dd))
      }
      res[1, ] <- rr2$parameters$mean
      res[2, ] <- rr2$parameters$variance$sigmasq
      res[3, ] <- rr2$parameters$pro
      
      ## bimodality index (BI)
      smd <- abs(res[1, 2] - res[1, 1]) / sqrt((res[2, 2] + res[2, 1]) / 2)
      bi <- sqrt(res[3, 2] * (1 - res[3, 2])) * smd
      
      #classification
      mystatus[cc.ix] <- as.numeric(rr2$classification == 2)
      mystatus.proba[cc.ix] <- rr2$z[ , 2, drop=TRUE]
      return(list("status"=mystatus, "status1.proba"=mystatus.proba, "gaussians"=res, "BIC"=mybic,  "BI"=bi, "x"=dd))
    } else { stop("Not enough data!") }
  }