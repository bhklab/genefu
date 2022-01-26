#' @title Function to fit a Single Sample Predictor (SSP) as in 
#'   Perou, Sorlie, Hu, and Parker publications
#'
#' @description
#' This function fits the Single Sample Predictor (SSP) as published
#'   in Sorlie et al 2003, Hu et al 2006 and Parker et al 2009. This model 
#'   is actually a nearest centroid classifier where the centroids
#'   representing the breast cancer molecular subtypes are identified
#'   through hierarchical clustering using an "intrinsic gene list".
#'
#' @usage
#' intrinsic.cluster(data, annot, do.mapping = FALSE, mapping,
#'   std = c("none", "scale", "robust"), rescale.q = 0.05, intrinsicg,
#'   number.cluster = 3, mins = 5, method.cor = c("spearman", "pearson"),
#'   method.centroids = c("mean", "median", "tukey"), filen, verbose = FALSE)
#'
#' @param data Matrix of gene expressions with samples in rows and probes in
#'   columns, dimnames being properly defined.
#' @param annot    Matrix of annotations with at least one column named 
#'   "EntrezGene.ID", dimnames being properly defined.
#' @param do.mapping    TRUE if the mapping through Entrez Gene ids must 
#'   be performed (in case of ambiguities, the most variant probe is kept 
#'   for each gene), FALSE otherwise.
#' @param mapping Matrix with columns "EntrezGene.ID" and "probe" used to 
#'   force the mapping such that the probes are not selected based on their 
#'   variance.
#' @param std Standardization of gene expressions: scale for traditional 
#'   standardization based on mean and standard deviation, robust for 
#'   standardization based on the 0.025 and 0.975 quantiles, none to keep 
#'   gene expressions unchanged.
#' @param rescale.q Proportion of expected outliers for (robust) rescaling 
#'   the gene expressions.
#' @param intrinsicg Intrinsic gene lists. May be specified by the user as 
#'   a matrix with at least 2 columns named probe and EntrezGene.ID for the 
#'   probe names and the corresponding Entrez Gene ids. The intrinsic gene lists 
#'   published by Sorlie et al. 2003, Hu et al. 2006 and Parker et al. 2009 are 
#'   stored in ssp2003, ssp2006 and pam50 respectively.
#' @param number.cluster The number of main clusters to be identified by 
#'   hierarchical clustering.
#' @param mins The minimum number of samples to be in a main cluster.
#' @param method.cor Correlation coefficient used to identified the nearest
#'   centroid. May be spearman or pearson.
#' @param method.centroids Method to compute a centroid from gene expressions 
#'   of a cluster of samples: mean, median or tukey (Tukey's Biweight Robust Mean).
#' @param filen Name of the csv file where the subtype clustering model must be 
#'   stored.
#' @param verbose TRUE to print informative messages, FALSE otherwise.
#'
#' @return
#' A list with items:
#' - model: Single Sample Predictor
#' - subtype: Subtypes identified by the SSP. For published intrinsic gene lists,
#'   subtypes can be either "Basal", "Her2", "LumA", "LumB" or "Normal".
#' - subtype.proba: Probabilities to belong to each subtype estimated from the
#'   correlations to each centroid.
#' - cor: Correlation coefficient to each centroid.
#'
#' @references
#' T. Sorlie and R. Tibshirani and J. Parker and T. Hastie and J. S. Marron and 
#'   A. Nobel and S. Deng and H. Johnsen and R. Pesich and S. Geister and J. 
#'   Demeter and C. Perou and P. E. Lonning and P. O. Brown and A. L. Borresen-Dale 
#'   and D. Botstein (2003) "Repeated Observation of Breast Tumor Subtypes in 
#'   Independent Gene Expression Data Sets", Proceedings of the National Academy of 
#'   Sciences, 1(14):8418-8423
#' 
#' Hu, Zhiyuan and Fan, Cheng and Oh, Daniel and Marron, JS and He, Xiaping and 
#'   Qaqish, Bahjat and Livasy, Chad and Carey, Lisa and Reynolds, Evangeline and 
#'   ressler, Lynn and Nobel, Andrew and Parker, Joel and Ewend, Matthew and Sawyer,
#'   Lynda and Wu, Junyuan and Liu, Yudong and Nanda, Rita and Tretiakova, Maria 
#'   and Orrico, Alejandra and Dreher, Donna and Palazzo, Juan and Perreard, 
#'   Laurent and Nelson, Edward and Mone, Mary and Hansen, Heidi and Mullins, 
#'   Michael and Quackenbush, John and Ellis, Matthew and Olopade, Olufunmilayo and 
#'   Bernard, Philip and Perou, Charles (2006) "The molecular portraits of breast 
#'   tumors are conserved across microarray platforms", BMC Genomics, 7(96)
#'
#' Parker, Joel S. and Mullins, Michael and Cheang, Maggie C.U. and Leung, 
#'   Samuel and Voduc, David and Vickery, Tammi and Davies, Sherri and Fauron, 
#'   Christiane and He, Xiaping and Hu, Zhiyuan and Quackenbush, John F. and 
#'   Stijleman, Inge J. and Palazzo, Juan and Marron, J.S. and Nobel, Andrew B. 
#'   and Mardis, Elaine and Nielsen, Torsten O. and Ellis, Matthew J. and Perou, 
#'   Charles M. and Bernard, Philip S. (2009) "Supervised Risk Predictor of Breast 
#'   Cancer Based on Intrinsic Subtypes", Journal of Clinical Oncology, 
#'   27(8):1160-1167
#'
#' @seealso
#' [genefu::subtype.cluster], [genefu::intrinsic.cluster.predict], [genefu::ssp2003], [genefu::ssp2006], [genefu::pam50]
#'
#' @examples
#' # load SSP signature published in Sorlie et al. 2003
#' data(ssp2003)
#' # load NKI data
#' data(nkis)
#' # load VDX data
#' data(vdxs)
#' ssp2003.nkis <- intrinsic.cluster(data=data.nkis, annot=annot.nkis,
#'   do.mapping=TRUE, std="robust",
#'   intrinsicg=ssp2003$centroids.map[ ,c("probe", "EntrezGene.ID")],
#'   number.cluster=5, mins=5, method.cor="spearman",
#'   method.centroids="mean", verbose=TRUE)
#' str(ssp2003.nkis, max.level=1)
#'
#' @md
#' @export
intrinsic.cluster <- function(data, annot, do.mapping=FALSE,
        mapping, std=c("none", "scale", "robust"), rescale.q=0.05,
        intrinsicg, number.cluster=3, mins=5, 
        method.cor=c("spearman", "pearson"), 
        method.centroids=c("mean", "median", "tukey"), filen, verbose=FALSE) {

    if(missing(data) || missing(annot) || missing(intrinsicg)) { 
        stop("data, annot, and intrinsicg parameters must be specified") 
    }
    std <- match.arg(std)
    method.cor <- match.arg(method.cor)
    method.centroids <- match.arg(method.centroids)
    if (!is.matrix(data)) {
        data <- as.matrix(data)
    }

    ## mapping
    if (do.mapping) {
        ## mapping through EntrezGene.ID
        if (is.matrix(intrinsicg) || is.data.frame(intrinsicg)) {
            ## matrix of annotations instead of a list of EntrezGene.IDs
            tt <- as.character(intrinsicg[, "EntrezGene.ID"])
            names(tt) <- as.character(intrinsicg[, "probe"])
            intrinsicg <- tt
        }
        gt <- length(intrinsicg)
        ## EntrezGene.IDs should be numeric or character that can be tranformed into numeric
        ## remove (arbitrarily) duplicated gene ids
        intrinsicg <- intrinsicg[!duplicated(intrinsicg) & !is.na(intrinsicg)]
        gid <- as.character(annot[ , "EntrezGene.ID"])
        names(gid) <- as.character(annot[ , "probe"])
        gid.intrinsic <- as.character(intrinsicg)
        names(gid.intrinsic) <- paste("geneid", gid.intrinsic, sep=".")
        if(missing(mapping)) { # select the most variant probes using annotations
            # if multiple probes for one gene, keep the most variant
            rr <- geneid.map(geneid1=gid, data1=data, geneid2=gid.intrinsic, verbose=FALSE)
            nn <- match(rr$geneid2, gid.intrinsic)
            nn <- nn[!is.na(nn)]
            intrinsicg <- intrinsicg[nn]
            data <- rr$data1
        } else { # use a predefined mapping
            nn <- as.character(mapping[, "EntrezGene.ID"])
            # keep only intrinsic genes with mapped probes
            myx <- is.element(gid.intrinsic, nn)
            gid.intrinsic <- gid.intrinsic[myx]
            intrinsicg <- intrinsicg[myx, ]
            pp <- as.character(mapping[match(gid.intrinsic, nn), "probe"])
            myx <- is.element(pp, dimnames(data)[[2]])
            intrinsicg <- intrinsicg[myx, ]
            pp <- pp[myx]
            data <- data[, pp, drop=FALSE]
        }
    }
    else {
        if (is.matrix(intrinsicg) || is.data.frame(intrinsicg)) {
            ## matrix of annotations instead of a list of EntrezGene.IDs
            tt <- as.character(intrinsicg[ , "probe"])
            names(tt) <- as.character(intrinsicg[ , "probe"])
            intrinsicg <- tt
        }
        gt <- length(intrinsicg)
        if(all(!is.element(dimnames(data)[[2]], intrinsicg))) { stop("no probe in common -> annot or mapping parameters are necessary for the mapping process!") }
        ## no mapping are necessary
        intrinsicg <- intersect(intrinsicg, dimnames(data)[[2]])
        data <- data[, intrinsicg]
    }
    gm <- length(intrinsicg)
    if (gm == 0 || (sum(is.na(data)) / length(data)) > 0.9) {
        stop("none of the instrinsic genes are present or too many missing values!")
    }
    if (!is.null(names(intrinsicg))) {
        centroids.map <- cbind("probe"=dimnames(data)[[2]],
            "probe.centroids"=names(intrinsicg),
            "EntrezGene.ID"=as.character(annot[dimnames(data)[[2]],
            "EntrezGene.ID"]))
    } else {
        centroids.map <- cbind("probe"=dimnames(data)[[2]],
            "probe.centroids"=NA,
            "EntrezGene.ID"=as.character(annot[dimnames(data)[[2]],
            "EntrezGene.ID"]))
    }
    dimnames(centroids.map)[[1]] <- dimnames(data)[[2]]

    if (verbose) {
        message(sprintf("%i/%i probes are used for clustering", gm, gt))
    }

    switch(std,
    "scale"={
        data <- scale(data, center=TRUE, scale=TRUE)
        if (verbose) message("standardization of the gene expressions")
    }, 
    "robust"={
        data <- apply(data, 2, function(x) {
            (rescale(x, q=rescale.q, na.rm=TRUE) - 0.5) * 2
        })
        if(verbose) message("robust standardization of the gene expressions")
    },
    "none"={
        if(verbose) message("no standardization of the gene expressions")
    })

    ## compute the clustering and cut the dendrogram
    ## hierarchical clustering with correlation-based distance and average linkage
    hcl <- amap::hcluster(x=data, method="correlation", link="average")
    mins.ok <- FALSE
    nbc <- number.cluster
    while (!mins.ok) { ## until each cluster contains at least mins samples
        cl <- cutree(tree=hcl, k=nbc)
        tt <- table(cl)
        if (sum(tt >= mins) >= number.cluster) {
            if (nbc > number.cluster) { ## put NA for clusters with less than mins samples
                td <- names(tt)[tt < mins]
                cl[is.element(cl, td)] <- NA
                ## rename the clusters
                ucl <- sort(unique(cl))
                ucl <- ucl[!is.na(ucl)]
                cl2 <- cl
                for (i in 1:number.cluster) { 
                    cl2[cl == ucl[i] & !is.na(cl)] <- i 
                }
                cl <- cl2
            }
            mins.ok <- TRUE
        } else {
            nbc <- nbc + 1
            if (nbc > (nrow(data) - (number.cluster * mins))) {
                stop("clusters are too small (size < mins)!")
            }
        }
    }

    ## compute the centroids
    ## take the core samples in each cluster to compute the centroid
    ## not feasible due to low intra correlation within clusters!!!
    ## minimal pairwise cor of approx 0.3
    cl.centroids <- matrix(NA, nrow=ncol(data), ncol=number.cluster,
        dimnames=list(dimnames(data)[[2]], paste("cluster", 1:number.cluster, sep=".")))
    for (i in 1:ncol(cl.centroids)) {
        switch(method.centroids,
        "mean"={ 
            cl.centroids[ ,i] <- apply(X=data[cl == i & !is.na(cl), ,drop=FALSE], 
                MARGIN=2, FUN=mean, na.rm=TRUE, trim=0.025)
        }, 
        "median"={ 
            cl.centroids[, i] <- apply(X=data[cl == i & !is.na(cl), ,drop=FALSE], 
                MARGIN=2, FUN=median, na.rm=TRUE)
        }, 
        "tukey"={ cl.centroids[ ,i] <- apply(X=data[cl == i & !is.na(cl), ,drop=FALSE], 
            MARGIN=2, FUN=tbrm, na.rm=TRUE, C=9)
        })
    }

    ## apply the nearest centroid classifier to classify the samples again
    ncor <- t(apply(X=data, MARGIN=1, FUN=function(x, y, z) {
        return(cor(x, y, method=z, use="complete.obs"))
    }, y=cl.centroids, z=method.cor))
    ## negative correlations are truncated to zero since they have no meaning for subtypes identification
    nproba <- t(apply(X=ncor, MARGIN=1, FUN=function(x) {
        x[x < 0] <- 0; return(x / sum(x, na.rm=TRUE)); 
    }))
    dimnames(ncor) <- dimnames(nproba) <- list(dimnames(data)[[1]], dimnames(cl.centroids)[[2]])
    ncl <- apply(X=ncor, MARGIN=1, FUN=function(x) {
        order(x, decreasing=TRUE, na.last=TRUE)[1]
    })
    ncl <- dimnames(cl.centroids)[[2]][ncl]
    names(ncl) <- dimnames(data)[[1]]

    if (!missing(filen)) {
        #save model parameters in a csv file for reuse
        write(x=sprintf("# Benjamin Haibe-Kains. All rights reserved."), sep="", file=paste(filen, "csv", sep="."))
        write(x=sprintf("# method.cor: %s", method.cor), sep="", append=TRUE, file=paste(filen, "csv", sep="."))
        write(x=sprintf("# method.centroids: %s", method.centroids), sep="", append=TRUE, file=paste(filen, "csv", sep="."))
        write(x=sprintf("# std: %s", std), sep="", append=TRUE, file=paste(filen, "csv", sep="."))
        write(x=sprintf("# rescale.q: %s", rescale.q), sep="", append=TRUE, file=paste(filen, "csv", sep="."))
        write(x=sprintf("# mins: %i", mins), sep="", append=TRUE, file=paste(filen, "csv", sep="."))
        ## centroids
        mycent <- t(cl.centroids)
        for(i in 1:nrow(mycent)) { 
            write(x=sprintf("# centroid.%s: %s", dimnames(mycent)[[1]][i], 
                paste(mycent[i, ], collapse=" ")), sep="", append=TRUE, 
                file=paste(filen, "csv", sep=".")) 
        }
        ## centroids.map
        write(paste("\"", dimnames(centroids.map)[[2]], "\"", collapse=",", sep=""), sep="", append=TRUE, file=paste(filen, "csv", sep="."))
        write.table(centroids.map, sep=",", col.names=FALSE, row.names=FALSE, file=paste(filen, "csv", sep="."), append=TRUE)
    }

    return(list(
        "model"=list("method.cor"=method.cor, 
            "method.centroids"=method.centroids, "std"=std,
            "rescale.q"=rescale.q,  "mins"=mins, "centroids"=cl.centroids,
            "centroids.map"=centroids.map),
        "subtype"=ncl,
        "subtype.proba"=nproba,
        "cor"=ncor
        )
    )
}