if(getRversion() >= "2.15.1")  utils::globalVariables("sigOvcTCGA")

`ovcTCGA` <-
function(data, annot, gmap=c("entrezgene", "ensembl_gene_id", "hgnc_symbol", "unigene"), do.mapping=FALSE, verbose=FALSE) {
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