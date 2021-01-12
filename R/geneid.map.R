#' @title Function to find the common genes between two datasets or a dataset and
#'   a gene list
#'
#' @description
#' This function allows for fast mapping between two datasets or a dataset and a gene
#'   list. The mapping process is performed using Entrez Gene id as reference. In case of
#'   ambiguities (several probes representing the same gene), the most variant probe is
#'   selected.
#'
#' @usage
#' geneid.map(geneid1, data1, geneid2, data2, verbose = FALSE)
#'
#' @param geneid1 First vector of Entrez Gene ids. The name of the vector cells must
#'   be the name of the probes in the dataset data1.
#' @param data1	First dataset with samples in rows and probes in columns. The dimnames
#'   must be properly defined.
#' @param geneid2 Second vector of Entrez Gene ids. The name of the vector cells must
#'   be the name of the probes in the dataset data1 if it is not missing, proper names must be assigned otherwise.
#' @param data2	First dataset with samples in rows and probes in columns. The dimnames
#'   must be properly defined. It may be missing.
#' @param verbose TRUE to print informative messages, FALSE otherwise.
#'
#'
#' @return
#' A list with items:
#' - geneid1 Mapped gene list from geneid1.
#' - data1 Mapped dataset from data1.
#' - geneid2 Mapped gene list from geneid2.
#' - data2 Mapped dataset from data2.
#'
#' @note
#' It is mandatory that the names of geneid1 and geneid2 must be the probe names
#'   of the microarray platform.
#'
#' @examples
#' # load NKI data
#' data(nkis)
#' nkis.gid <- annot.nkis[ ,"EntrezGene.ID"]
#' names(nkis.gid) <- dimnames(annot.nkis)[[1]]
#' # load GGI signature
#' data(sig.ggi)
#' ggi.gid <- sig.ggi[ ,"EntrezGene.ID"]
#' names(ggi.gid) <- as.character(sig.ggi[ ,"probe"])
#' # mapping through Entrez Gene ids of NKI and GGI signature
#' res <- geneid.map(geneid1=nkis.gid, data1=data.nkis,
#'   geneid2=ggi.gid, verbose=FALSE)
#' str(res)
#'
#' @md
#' @export
geneid.map <-
function(geneid1, data1, geneid2, data2, verbose=FALSE) {

	nn <- names(geneid1)
	geneid1 <- as.character(geneid1)
	names(geneid1) <- nn
	nn <- names(geneid2)
	geneid2 <- as.character(geneid2)
	names(geneid2) <- nn
	if(is.null(names(geneid1))) { names(geneid1) <- dimnames(data1)[[2]] }
	if(!missing(data2) && is.null(names(geneid2))) { names(geneid2) <- dimnames(data2)[[2]] }
	if(!missing(data1) && !missing(geneid1) && !missing(geneid2)) {
		## remove probes without any measurements
		na.ix <- apply(data1, 2, function(x) { return(all(is.na(x))) })
		data1 <- data1[ , !na.ix, drop=FALSE]
		geneid1 <- geneid1[!na.ix]
	} else { stop("data1, geneid1 and geneid2 parameters are mandatory!") }
	if(!missing(data2)) {
		## remove probes without any measurements
		na.ix <- apply(data2, 2, function(x) { return(all(is.na(x))) })
		data2 <- data2[ , !na.ix, drop=FALSE]
		geneid2 <- geneid2[!na.ix]
	} else { data2 <- NULL }

	gix1 <- !is.na(geneid1)
	gix2 <- !is.na(geneid2)

	geneid.common <- intersect(geneid1[gix1], geneid2[gix2])
	if(length(geneid.common) == 0) {
		warning("no gene ids in common!")
		return(list("geneid1"=NA, "data1"=NA, "geneid2"=NA, "data2"=NA))
	}

	## dataset1
	## probes corresponding to common gene ids
	gg <- names(geneid1)[is.element(geneid1, geneid.common)]
	gid <- geneid1[is.element(geneid1, geneid.common)]
	## duplicated gene ids
	gid.dupl <- unique(gid[duplicated(gid)])
	gg.dupl <- names(geneid1)[is.element(geneid1, gid.dupl)]
	## unique gene ids
	gid.uniq <- gid[!is.element(gid, gid.dupl)]
	gg.uniq <- names(geneid1)[is.element(geneid1, gid.uniq)]
	## data corresponding to unique gene ids
	datat <- data1[ ,gg.uniq,drop=FALSE]
	## data for duplicated gene ids
	if(length(gid.dupl) > 0) {
		if(verbose) { message("\ndataset1 duplicates...") }
		## compute the standard deviation with a penalization on the number of missing values
		## this should avoid selecting the most variant probe with a lot of missing values
		pena <- apply(X=data1[ , gg.dupl, drop=FALSE], MARGIN=2, FUN=function(x) { return(sum(is.na(x))) })
		pena <- log((nrow(data1) + 1) / (pena + 1)) + 1
		#pena <- 1
		sdr <- drop(apply(X=data1[ , gg.dupl, drop=FALSE], MARGIN=2, FUN=sd, na.rm=TRUE)) * pena
		mysd <- cbind("probe"=gg.dupl, "gid"=geneid1[gg.dupl], "sd"=sdr)
		mysd <- mysd[order(as.numeric(mysd[ , "sd"]), decreasing=TRUE, na.last=TRUE), , drop=FALSE]
		mysd <- mysd[!duplicated(mysd[ , "gid"]), , drop=FALSE]
		datat <- cbind(datat, data1[ , mysd[ , "probe"], drop=FALSE])
	}
	data1 <- datat
	geneid1 <- geneid1[dimnames(data1)[[2]]]

	#dataset2
	if(is.null(data2)) {
		#keep arbitrarily the first occurence of each duplicated geneid
		geneid2 <- geneid2[!duplicated(geneid2) & is.element(geneid2, geneid.common)]
	}
	else {
		## probes corresponding to common gene ids
		gg <- names(geneid2)[is.element(geneid2, geneid.common)]
		gid <- geneid2[is.element(geneid2, geneid.common)]
		## duplicated gene ids
		gid.dupl <- unique(gid[duplicated(gid)])
		gg.dupl <- names(geneid2)[is.element(geneid2, gid.dupl)]
		## unique gene ids
		gid.uniq <- gid[!is.element(gid, gid.dupl)]
		gg.uniq <- names(geneid2)[is.element(geneid2, gid.uniq)]
		## data corresponding to unique gene ids
		datat <- data2[ ,gg.uniq,drop=FALSE]
		## data for duplicated gene ids
		if(length(gid.dupl) > 0) {
			if(verbose) { message("\ndataset2 duplicates...") }
			## compute the standard deviation with a penalization on the number of missing values
			## this should avoid selecting the most variant probe with a lotof missing values
			pena <- apply(X=data2[ , gg.dupl, drop=FALSE], MARGIN=2, FUN=function(x) { return(sum(is.na(x))) })
			pena <- log((nrow(data2) + 1) / (pena + 1)) + 1
			#pena <- 1
			sdr <- drop(apply(X=data2[ , gg.dupl, drop=FALSE], MARGIN=2, FUN=sd, na.rm=TRUE)) * pena
			mysd <- cbind("probe"=gg.dupl, "gid"=geneid2[gg.dupl], "sd"=sdr)
			mysd <- mysd[order(as.numeric(mysd[ , "sd"]), decreasing=TRUE, na.last=TRUE), , drop=FALSE]
			mysd <- mysd[!duplicated(mysd[ , "gid"]), , drop=FALSE]
			datat <- cbind(datat, data2[ , mysd[ , "probe"], drop=FALSE])
		}
		data2 <- datat
		geneid2 <- geneid2[dimnames(data2)[[2]]]
	}

	#same order for the two datasets
	rix <- match(geneid2, geneid1)
	geneid1 <- geneid1[rix]
	data1 <- data1[ ,rix,drop=FALSE]
	return(list("geneid1"=geneid1, "data1"=data1, "geneid2"=geneid2, "data2"=data2))
}