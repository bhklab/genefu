`ovcAngiogenic` <-
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
        mymapping <- c("mapped"=nrow(sigt), "total"=nrow(sigOvcTCGA))
        myprobe <- data.frame("probe"=pold, "gene.map"=annot[ ,gmap], "new.probe"=pold2)
    } else {
        gix <- intersect(rownames(sigOvcAngiogenic), colnames(data))
        if(length(gix) < 2) { stop("data do not contain enough gene from the ovcTCGA signature!") }
        data <- data[ ,gix,drop=FALSE]
        annot <- annot[gix, ,drop=FALSE]
        mymapping <- c("mapped"=length(gix), "total"=nrow(sigOvcTCGA))
        myprobe <- data.frame("probe"=gix, "gene.map"=annot[ ,gmap], "new.probe"=gix)
        sigt <- sigOvcAngiogenic[gix, ,drop=FALSE]
    }
    
    ss <- genefu::sig.score(x=data.frame("probe"=colnames(data), "EntrezGene.ID"=annot[ ,gmap], "coefficient"=sigt[ ,"weight"]), data=data, annot=annot, do.mapping=FALSE, signed=TRUE)$score
    ## rescale only with the high grade, late stage, serous (hgs) patients
    rr <- genefu::rescale(ss[hgs], q=0.05, na.rm=TRUE)
    ## rescale the whole dataset
    ## rescale the whole dataset
    pscore <- ((ss - attributes(rr)$q1) / (attributes(rr)$q2 - attributes(rr)$q1) - 0.5) * 2
    emclust.ts <- mclust::estep(modelName="E", data=pscore, parameters=modelOvcAngioganic)
    dimnames(emclust.ts$z) <- list(names(pscore), c("Angiogenic.proba", "nonAngiogenic.proba"))
    class.ts <- mclust::map(emclust.ts$z, warn=FALSE)
    names(class.ts) <- names(pscore)
    sbt.ts <- class.ts
    sbt.ts[class.ts == 1] <- "Angiogenic"
    sbt.ts[class.ts == 2] <- "nonAngiogenic"
    sbts <- data.frame("subtype. Angiogenic.score"=pscore, "subtype.Angiogenic"=sbt.ts, emclust.ts$z)
    prisk <- as.numeric(pscore <= 0.5)
	names(prisk) <- names(pscore) <- rownames(data)
	return (list("score"=pscore, "risk"=prisk, "mapping"=mymapping, "probe"=myprobe, "subtype"=sbts))
}