`.ovcSigs` <-
function(sigs=c("bentink2012_angiogenic", "crijns2009_sig", "yoshihara2010_sig", "spentzos2011_sig", "tcga2011_sig")) {
    for(i in 1:length(sigs)) {
        sig <- read.csv(system.file(file.path("extdata", sprintf("%s.csv", sigs[i])), package="genefu"), stringsAsFactors=FALSE)
        ## annotations
        ensembl.db <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
        switch(sigs[i],
        "bentink2012_angiogenic"={
            ss <- "illumina_humanwg_6_v2"
            gid <- as.character(sig[ ,"Probe_Id"])
            gene.an <- biomaRt::getBM(attributes=c(ss, "entrezgene", "ensembl_gene_id", "hgnc_symbol", "unigene", "description", "chromosome_name", "start_position", "end_position", "strand", "band"), filters=ss, values=sort(unique(gid)), mart=ensembl.db)
            gene.an[gene.an == "" | gene.an == " "] <- NA
            gene.an <- gene.an[!is.na(gene.an[ , ss]) & !duplicated(gene.an[ , ss]) & is.element(gene.an[ , ss], gid), , drop=FALSE]
            annot <- data.frame(matrix(NA, nrow=nrow(sig), ncol=ncol(gene.an), dimnames=list(gid, colnames(gene.an))))
            annot[match(gene.an[ , ss], gid), colnames(gene.an)] <- gene.an
            annot <- data.frame("probe"=gid, annot, "weight"=as.numeric(sig[ ,"weights"]))
            sigOvcAngiogenic <- annot
            #save(list="sigAngiogenic", compress=TRUE, file=file.path(system.file(package="genefu"), "data", "sigAngiogenic.rda"))
            save(list="sigOvcAngiogenic", compress=TRUE, file="sigOvcAngiogenic.rda")
        },
        "crijns2009_sig"={
            ss <- "hgnc_symbol"
            gid <- as.character(sig[ ,"Gene.Id"])
            gene.an <- biomaRt::getBM(attributes=c(ss, "entrezgene", "ensembl_gene_id", "unigene", "description", "chromosome_name", "start_position", "end_position", "strand", "band"), filters=ss, values=sort(unique(gid)), mart=ensembl.db)
            gene.an[gene.an == "" | gene.an == " "] <- NA
            gene.an <- gene.an[!is.na(gene.an[ , ss]) & !duplicated(gene.an[ , ss]) & is.element(gene.an[ , ss], gid), , drop=FALSE]
            annot <- data.frame(matrix(NA, nrow=nrow(sig), ncol=ncol(gene.an), dimnames=list(gid, colnames(gene.an))))
            annot[match(gene.an[ , ss], gid), colnames(gene.an)] <- gene.an
            annot <- data.frame("probe"=gid, annot, "weight"=as.numeric(sig[ ,"Weights"]))
            sigOvcCrijns <- annot
            save(list="sigOvcCrijns", compress=TRUE, file="sigOvcCrijns.rda")
        },
        "yoshihara2010_sig"={
            ss <- "refseq_mrna"
            gid <- as.character(sig[ ,"GenBank"])
            gene.an <- biomaRt::getBM(attributes=c(ss, "entrezgene", "ensembl_gene_id", "unigene", "description", "chromosome_name", "start_position", "end_position", "strand", "band"), filters=ss, values=sort(unique(gid)), mart=ensembl.db)
            gene.an[gene.an == "" | gene.an == " "] <- NA
            gene.an <- gene.an[!is.na(gene.an[ , ss]) & !duplicated(gene.an[ , ss]) & is.element(gene.an[ , ss], gid), , drop=FALSE]
            annot <- data.frame(matrix(NA, nrow=nrow(sig), ncol=ncol(gene.an), dimnames=list(gid, colnames(gene.an))))
            annot[match(gene.an[ , ss], gid), colnames(gene.an)] <- gene.an
            gene.an <- biomaRt::getBM(attributes=c(ss, "hgnc_symbol"), filters=ss, values=sort(unique(gid)), mart=ensembl.db)
            gene.an[gene.an == "" | gene.an == " "] <- NA
            gene.an <- gene.an[!is.na(gene.an[ , ss]) & !duplicated(gene.an[ , ss]) & is.element(gene.an[ , ss], gid), , drop=FALSE]
            annot2 <- data.frame(matrix(NA, nrow=nrow(sig), ncol=ncol(gene.an), dimnames=list(gid, colnames(gene.an))))
            annot2[match(gene.an[ , ss], gid), colnames(gene.an)] <- gene.an
            annot <- data.frame("probe"=gid, annot2[ ,-1,drop=FALSE], annot, "weight"=as.numeric(sig[ ,"beta_ridge"]))
            sigOvcYoshihara <- annot
            save(list="sigOvcYoshihara", compress=TRUE, file="sigOvcYoshihara.rda")
        },
        "spentzos2011_sig"={
            ss <- "affy_hg_u133a"
            gid <- as.character(sig[ ,"probeset"])
            gene.an <- biomaRt::getBM(attributes=c(ss, "entrezgene", "hgnc_symbol", "ensembl_gene_id", "unigene", "description", "chromosome_name", "start_position", "end_position", "strand", "band"), filters=ss, values=sort(unique(gid)), mart=ensembl.db)
            gene.an[gene.an == "" | gene.an == " "] <- NA
            gene.an <- gene.an[!is.na(gene.an[ , ss]) & !duplicated(gene.an[ , ss]) & is.element(gene.an[ , ss], gid), , drop=FALSE]
            annot <- data.frame(matrix(NA, nrow=nrow(sig), ncol=ncol(gene.an), dimnames=list(gid, colnames(gene.an))))
            annot[match(gene.an[ , ss], gid), colnames(gene.an)] <- gene.an
            annot <- data.frame("probe"=gid, annot, "weight"=sig[ ,"weight"])
            sigOvcSpentzos <- annot
            save(list="sigOvcSpentzos", compress=TRUE, file="sigOvcSpentzos.rda")
        },
        "tcga2011_sig"={
            ss <- "entrezgene"
            gid <- as.character(sig[ ,"Entrez.Id"])
            gene.an <- biomaRt::getBM(attributes=c(ss, "hgnc_symbol", "ensembl_gene_id", "unigene", "description", "chromosome_name", "start_position", "end_position", "strand", "band"), filters=ss, values=sort(unique(gid)), mart=ensembl.db)
            gene.an[gene.an == "" | gene.an == " "] <- NA
            gene.an <- gene.an[!is.na(gene.an[ , ss]) & !duplicated(gene.an[ , ss]) & is.element(gene.an[ , ss], gid), , drop=FALSE]
            annot <- data.frame(matrix(NA, nrow=nrow(sig), ncol=ncol(gene.an), dimnames=list(paste("geneid", gid, sep="."), colnames(gene.an))))
            annot[match(gene.an[ , ss], gid), colnames(gene.an)] <- gene.an
            annot <- data.frame("probe"=gid, annot, sig[ ,c("Gene.set", "beta", "p.value")])
            sigOvcTCGA <- annot
            save(list="sigOvcTCGA", compress=TRUE, file="sigOvcTCGA.rda")
        })
    }
}