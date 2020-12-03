#' @title Function to compute the subtype scores and risk classifications for 
#'   the prognostic signature published by Yoshihara et al.
#'
#' @description
#' This function computes subtype scores and risk classifications from gene 
#'   expression values following the algorithm developed by Yoshihara et al, 
#'   for prognosis in ovarian cancer.
#'
#' @usage
#' ovcYoshihara(data, annot, hgs, 
#'   gmap = c("entrezgene", "ensembl_gene_id", "hgnc_symbol", "unigene", "refseq_mrna"), 
#'   do.mapping = FALSE, verbose = FALSE)
#'
#' @param data	Matrix of gene expressions with samples in rows and probes in 
#'   columns, dimnames being properly defined.
#' @param annot	Matrix of annotations with one column named as gmap, dimnames 
#'   being properly defined.
#' @param hgs vector of booleans with TRUE represents the ovarian cancer 
#'   patients who have a high grade, late stage, serous tumor, FALSE otherwise. 
#'   This is particularly important for properly rescaling the data. If hgs is 
#'   missing, all the patients will be used to rescale the subtype score.
#' @param gmap character string containing the biomaRt attribute to use for 
#'   mapping if do.mapping=TRUE
#' @param do.mapping TRUE if the mapping through Entrez Gene ids must be 
#'   performed (in case of ambiguities, the most variant probe is kept for 
#'   each gene), FALSE otherwise.
#' @param verbose TRUE to print informative messages, FALSE otherwise.
#'
#' @return
#' A list with items:
#' - score: Continuous signature scores.
#' - risk: Binary risk classification, 1 being high risk and 0 being low risk.
#' - mapping: Mapping used if necessary.
#' - probe: If mapping is performed, this matrix contains the correspondence 
#'   between the gene list (aka signature) and gene expression data.
#'
#' @references
#' Yoshihara K, Tajima A, Yahata T, Kodama S, Fujiwara H, Suzuki M, Onishi Y, 
#'   Hatae M, Sueyoshi K, Fujiwara H, Kudo, Yoshiki, Kotera K, Masuzaki H, 
#'   Tashiro H, Katabuchi H, Inoue I, Tanaka K (2010) "Gene expression profile 
#'   for predicting survival in advanced-stage serous ovarian cancer across two 
#'   independent datasets", PloS one, 5(3):e9615.
#'
#' @seealso
#' [genefu::sigOvcYoshihara]
#' 
#' @examples
#' # load the ovcYoshihara signature
#' data(sigOvcYoshihara)
#' # load NKI dataset
#' data(nkis)
#' colnames(annot.nkis)[is.element(colnames(annot.nkis), "EntrezGene.ID")] <- "entrezgene"
#' # compute relapse score
#' ovcYoshihara.nkis <- ovcYoshihara(data=data.nkis, 
#'   annot=annot.nkis, gmap="entrezgene", do.mapping=TRUE)
#' table(ovcYoshihara.nkis$risk)
#'
#' @md
#' @export
if(getRversion() >= "2.15.1")  utils::globalVariables("sigOvcYoshihara")

ovcYoshihara <-
function(data, annot, hgs, gmap=c("entrezgene", "ensembl_gene_id", "hgnc_symbol", "unigene", "refseq_mrna"), do.mapping=FALSE, verbose=FALSE) {
    gmap <- match.arg(gmap)
    if(missing(hgs)) { hgs <- rep(TRUE, nrow(data)) }
    if(do.mapping) {
        if(!is.element(gmap, colnames(annot))) { stop("gmap is not a column of annot!") }
        if(verbose) { message("the most variant probe is selected for each gene") }
        sigt <- sigOvcYoshihara[order(abs(sigOvcYoshihara[ ,"weight"]), decreasing=FALSE), ,drop=FALSE]
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
        mymapping <- c("mapped"=nrow(sigt), "total"=nrow(sigOvcYoshihara))
        myprobe <- data.frame("probe"=pold, "gene.map"=annot[ ,gmap], "new.probe"=pold2)
    } else {
        gix <- intersect(rownames(sigOvcYoshihara), colnames(data))
        if(length(gix) < 2) { stop("data do not contain enough gene from the ovcTCGA signature!") }
        data <- data[ ,gix,drop=FALSE]
        annot <- annot[gix, ,drop=FALSE]
        mymapping <- c("mapped"=length(gix), "total"=nrow(sigOvcYoshihara))
        myprobe <- data.frame("probe"=gix, "gene.map"=annot[ ,gmap], "new.probe"=gix)
        sigt <- sigOvcYoshihara[gix, ,drop=FALSE]
    }
    ## transform the gene expression in Z-scores
    data <- scale(data)
    pscore <- genefu::sig.score(x=data.frame("probe"=colnames(data), "EntrezGene.ID"=annot[ ,gmap], "coefficient"=sigt[ ,"weight"]), data=data, annot=annot, do.mapping=FALSE, signed=FALSE)$score
    prisk <- as.numeric(pscore > median(pscore, na.rm=TRUE))
	names(prisk) <- names(pscore) <- rownames(data)
	return (list("score"=pscore, "risk"=prisk, "mapping"=mymapping, "probe"=myprobe))
}