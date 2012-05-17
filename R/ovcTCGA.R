`pik3cags` <-
function(data, annot, gmap=c("entrezgene", "ensembl_gene_id", "hgnc_symbol", "unigene"), do.mapping=FALSE, verbose=FALSE) {

    gmap <- match.arg(gmap)
    if(do.mapping) {
        if(!is.element(gmap, colnames(annot))) { stop("gmap is not a column of annot!") }
        if(verbose) { message("the most variant probe is selected for each gene") }
        sigOvcTCGA <- sigOvcTCGA[order(sigOvcTCGA[ ,"p.value"], decreasing=TRUE), ,drop=FALSE]
        sigt <- sigOvcTCGA[!duplicated(sigOvcTCGA[ ,gmap]), ,drop=FALSE]
        gid2 <- sigt[ ,gmap]
        names(gid2) <- rownames(sigt)
        gid1 <- annot[ ,gmap]
        names(gid1) <- colnames(data)
        rr <- geneid.map(geneid1=gid1, data1=data, geneid2=gid2)
        data <- rr$data1
        annot <- annot[colnames(data), ,drop=FALSE]
        colnames(data) <- rownames(annot) <- rownames(sigt) <- paste("geneid", annot[ ,gmap], sep=".")
    } else {
        gix <- intersect(rownames(sigOvcTCGA), colnames(data))
        if(length(gix) < 2) { stop("data do not contain enough gene from the ovcTCGA signature!") }
        data <- data[ ,gix,drop=FALSE]
        annot <- annot[gix, ,drop=FALSE]
        mymapping <- c("mapped"=length(gix), "total"=nrow(sigOvcTCGA))
        myprobe <- data.frame("probe"=gix, "gene.map"=annot[ ,gmap], "new.probe"=gix)
        sigt <- sigOvcTCGA[gix, ,drop=FALSE]
    }
    myprobe <- data.frame("probe"=pold, )
    pscore <- apply(data, 1, function(x, y) { return(t.test(x ~ y)$statistic) }, y=as.numeric(sigt[ ,"beta"] < 0))
	prisk <- pscore >= 0
	names(prisk) <- names(pscore) <- rownames(data)
	return (list("score"=pscore, "risk"=prisk, "mapping"=, "probe"=))
}