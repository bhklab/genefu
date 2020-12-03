#' @title Function to compute the subtype scores and risk classifications 
#'   for the angiogenic molecular subtype in ovarian cancer
#'
#' @description
#' This function computes subtype scores and risk classifications from 
#'   gene expression values following the algorithm developed by Bentink, 
#'   Haibe-Kains et al. to identify the angiogenic molecular subtype in 
#'   ovarian cancer.
#'
#' @usage
#' ovcAngiogenic(data, annot, hgs, 
#' gmap = c("entrezgene", "ensembl_gene_id", "hgnc_symbol", "unigene"), 
#' do.mapping = FALSE, verbose = FALSE)
#'
#' @param data Matrix of gene expressions with samples in rows and probes in 
#'   columns, dimnames being properly defined.
#' @param annot Matrix of annotations with one column named as gmap, dimnames 
#'   being properly defined.
#' @param hgs vector of booleans with TRUE represents the ovarian cancer 
#'   patients who have a high grade, late stage, serous tumor, FALSE otherwise. This is particularly important for properly rescaling the data. If hgs is missing, all the patients will be used to rescale the subtype score.
#' @param gmap character string containing the biomaRt attribute to use for 
#'   mapping if do.mapping=TRUE
#' @param do.mapping TRUE if the mapping through Entrez Gene ids must be 
#'   performed (in case of ambiguities, the most variant probe is kept for each gene), FALSE otherwise.
#' @param verbose TRUE to print informative messages, FALSE otherwise.
#'
#'
#' @return
#' A list with items:
#' - score: Continuous signature scores.
#' - risk: Binary risk classification, 1 being high risk and 0 being low risk.
#' - mapping: Mapping used if necessary.
#' - probe: If mapping is performed, this matrix contains the correspondence 
#'   between the gene list (aka signature) and gene expression data.
#' - subtype: data frame reporting the subtype score, maximum likelihood 
#'   classification and corresponding subtype probabilities.
#'
#' @references
#' Bentink S, Haibe-Kains B, Risch T, Fan J-B, Hirsch MS, Holton K, Rubio R, 
#'   April C, Chen J, Wickham-Garcia E, Liu J, Culhane AC, Drapkin R, Quackenbush 
#'   JF, Matulonis UA (2012) "Angiogenic mRNA and microRNA Gene Expression 
#'   Signature Predicts a Novel Subtype of Serous Ovarian Cancer", PloS one, 
#'   7(2):e30269
#'
#' @seealso
#' [genefu::sigOvcAngiogenic]
#' 
#' @examples
#' # load the ovcAngiogenic signature
#' data(sigOvcAngiogenic)
#' # load NKI dataset
#' data(nkis)
#' colnames(annot.nkis)[is.element(colnames(annot.nkis), "EntrezGene.ID")] <- "entrezgene"
#' # compute relapse score
#' ovcAngiogenic.nkis <- ovcAngiogenic(data=data.nkis, annot=annot.nkis, 
#' gmap="entrezgene", do.mapping=TRUE)
#' table(ovcAngiogenic.nkis$risk)
#'
#' @md
#' @export
if(getRversion() >= "2.15.1")  utils::globalVariables(c("sigOvcAngiogenic","modelOvcAngiogenic"))

ovcAngiogenic <-
function(data, annot, hgs, gmap=c("entrezgene", "ensembl_gene_id", "hgnc_symbol", "unigene"), do.mapping=FALSE, verbose=FALSE) {
    gmap <- match.arg(gmap)
    if(missing(hgs)) { hgs <- rep(TRUE, nrow(data)) }
    if(do.mapping) {
        if(!is.element(gmap, colnames(annot))) { stop("gmap is not a column of annot!") }
        if(verbose) { message("the most variant probe is selected for each gene") }
        sigt <- sigOvcAngiogenic[order(abs(sigOvcAngiogenic[ ,"weight"]), decreasing=FALSE), ,drop=FALSE]
        sigt <- sigt[!duplicated(sigt[ ,gmap]), ,drop=FALSE]
        gid2 <- sigt[ ,gmap]
        names(gid2) <- rownames(sigt)
        gid1 <- annot[ ,gmap]
        names(gid1) <- colnames(data)
        rr <- geneid.map(geneid1=gid1, data1=data, geneid2=gid2)
        data <- rr$data1
        annot <- annot[colnames(data), ,drop=FALSE]
        sigt <- sigt[names(rr$geneid2), ,drop=FALSE]
        pold <- colnames(data)
        pold2 <- rownames(sigt)
        colnames(data) <- rownames(annot) <- rownames(sigt) <- paste("geneid", annot[ ,gmap], sep=".")
        mymapping <- c("mapped"=nrow(sigt), "total"=nrow(sigOvcAngiogenic))
        myprobe <- data.frame("probe"=pold, "gene.map"=annot[ ,gmap], "new.probe"=pold2)
    } else {
        gix <- intersect(rownames(sigOvcAngiogenic), colnames(data))
        if(length(gix) < 2) { stop("data do not contain enough gene from the ovcTCGA signature!") }
        data <- data[ ,gix,drop=FALSE]
        annot <- annot[gix, ,drop=FALSE]
        mymapping <- c("mapped"=length(gix), "total"=nrow(sigOvcAngiogenic))
        myprobe <- data.frame("probe"=gix, "gene.map"=annot[ ,gmap], "new.probe"=gix)
        sigt <- sigOvcAngiogenic[gix, ,drop=FALSE]
    }
    
    #data(modelOvcAngiogenic)
    ss <- genefu::sig.score(x=data.frame("probe"=colnames(data), "EntrezGene.ID"=annot[ ,gmap], "coefficient"=sigt[ ,"weight"]), data=data, annot=annot, do.mapping=FALSE, signed=TRUE)$score
    ## rescale only with the high grade, late stage, serous (hgs) patients
    rr <- genefu::rescale(ss[hgs], q=0.05, na.rm=TRUE)
    ## rescale the whole dataset
    pscore <- ((ss - attributes(rr)$q1) / (attributes(rr)$q2 - attributes(rr)$q1) - 0.5) * 2
    emclust.ts <- mclust::estep(modelName="E", data=pscore, parameters=modelOvcAngiogenic)
    dimnames(emclust.ts$z) <- list(names(pscore), c("Angiogenic.proba", "nonAngiogenic.proba"))
    class.ts <- mclust::map(emclust.ts$z, warn=FALSE)
    names(class.ts) <- names(pscore)
    sbt.ts <- class.ts
    sbt.ts[class.ts == 1] <- "Angiogenic"
    sbt.ts[class.ts == 2] <- "nonAngiogenic"
    sbts <- data.frame("subtype.score"=pscore, "subtype"=sbt.ts, emclust.ts$z)
    prisk <- as.numeric(sbts[ ,"subtype"] == "Angiogenic")
	names(prisk) <- names(pscore) <- rownames(data)
	return (list("score"=pscore, "risk"=prisk, "mapping"=mymapping, "probe"=myprobe, "subtype"=sbts))
}