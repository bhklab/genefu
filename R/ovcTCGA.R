#' @title Function to compute the prediction scores and risk classifications
#'   for the ovarian cancer TCGA signature
#'
#' @description
#' This function computes signature scores and risk classifications from gene
#'   expression values following the algorithm developed by the TCGA consortium
#'   for ovarian cancer.
#'
#' @usage
#' ovcTCGA(data, annot,
#'   gmap = c("entrezgene", "ensembl_gene_id", "hgnc_symbol", "unigene"),
#'   do.mapping = FALSE, verbose = FALSE)
#'
#' @param data	Matrix of gene expressions with samples in rows and probes in
#'   columns, dimnames being properly defined.
#' @param annot	Matrix of annotations with one column named as gmap, dimnames
#'   being properly defined.
#' @param gmap	character string containing the biomaRt attribute to use for
#'   mapping if do.mapping=TRUE
#' @param do.mapping TRUE if the mapping through Entrez Gene ids must be
#'   performed (in case of ambiguities, the most variant probe is kept for
#'   each gene), FALSE otherwise.
#' @param verbose TRUE to print informative messages, FALSE otherwise.
#'
#' @return
#' A list with items:
#' - score:	Continuous signature scores.
#' - risk: Binary risk classification, 1 being high risk and 0 being low risk.
#' - mapping: Mapping used if necessary.
#' - probe:	If mapping is performed, this matrix contains the correspondence
#'   between the gene list (aka signature) and gene expression data.
#'
#' @references
#' Bell D, Berchuck A, Birrer M et al. (2011) "Integrated genomic analyses of
#'   ovarian carcinoma", Nature, 474(7353):609-615
#'
#' @seealso
#' [genefu::sigOvcTCGA]
#'
#' @examples
#' # load the ovcTCGA signature
#' data(sigOvcTCGA)
#' # load NKI dataset
#' data(nkis)
#' colnames(annot.nkis)[is.element(colnames(annot.nkis), "EntrezGene.ID")] <- "entrezgene"
#' # compute relapse score
#' ovcTCGA.nkis <- ovcTCGA(data=data.nkis, annot=annot.nkis, gmap="entrezgene", do.mapping=TRUE)
#' table(ovcTCGA.nkis$risk)
#'
#' @md
#' @export
ovcTCGA <- function(data, annot, gmap=c("entrezgene", "ensembl_gene_id",
    "hgnc_symbol", "unigene"), do.mapping=FALSE, verbose=FALSE)
{
    gmap <- match.arg(gmap)
    if(do.mapping) {
        if(!is.element(gmap, colnames(annot))) { stop("gmap is not a column of annot!") }
        if(verbose) { message("the most variant probe is selected for each gene") }
        sigt <- sigOvcTCGA[order(sigOvcTCGA[ ,"p.value"], decreasing=TRUE), ,drop=FALSE]
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
        mymapping <- c("mapped"=nrow(sigt), "total"=nrow(sigOvcTCGA))
        myprobe <- data.frame("probe"=pold, "gene.map"=annot[ ,gmap], "new.probe"=pold2)
    } else {
        gix <- intersect(rownames(sigOvcTCGA), colnames(data))
        if(length(gix) < 2) { stop("data do not contain enough gene from the ovcTCGA signature!") }
        data <- data[ ,gix,drop=FALSE]
        annot <- annot[gix, ,drop=FALSE]
        mymapping <- c("mapped"=length(gix), "total"=nrow(sigOvcTCGA))
        myprobe <- data.frame("probe"=gix, "gene.map"=annot[ ,gmap], "new.probe"=gix)
        sigt <- sigOvcTCGA[gix, ,drop=FALSE]
    }
    pscore <- apply(data, 1, function(x, y) { return(t.test(x ~ y)$statistic) }, y=as.numeric(sigt[ ,"beta"] < 0))
	prisk <- as.numeric(pscore >= 0)
	names(prisk) <- names(pscore) <- rownames(data)
	return (list("score"=pscore, "risk"=prisk, "mapping"=mymapping, "probe"=myprobe))
}