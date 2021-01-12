#' @title Function to map a list of datasets through EntrezGene IDs in order to 
#'   get the union of the genes
#'
#' @description
#' This function maps a list of datasets through EntrezGene IDs in order to get 
#'   the union of the genes.
#'
#' @usage
#' map.datasets(datas, annots, do.mapping = FALSE, 
#'   mapping.coln = "EntrezGene.ID", mapping, verbose = FALSE)
#'
#' @param datas	List of matrices of gene expressions with samples in rows and 
#'   probes in columns, dimnames being properly defined.
#' @param annots	List of matrices of annotations with at least one column named 
#'   "EntrezGene.ID", dimnames being properly defined.
#' @param do.mapping	TRUE if the mapping through Entrez Gene ids must be 
#'   performed (in case of ambiguities, the most variant probe is kept for each 
#'   gene), FALSE otherwise.
#' @param mapping.coln	Name of the column containing the biological annotation 
#'   to be used to map the different datasets, default is "EntrezGene.ID".
#' @param mapping	Matrix with columns "EntrezGene.ID" and "probe.x" used to 
#'   force the mapping such that the probes of platform x are not selected based on 
#'   their variance.
#' @param verbose	TRUE to print informative messages, FALSE otherwise.
#'
#' @details
#' In case of several probes representing the same EntrezGene ID, the most 
#'   variant is selected if mapping is not specified. When a EntrezGene ID does not 
#'   exist in a specific dataset, NA values are introduced.
#'
#' @return
#' A list with items:
#' - datas: List of datasets (gene expression matrices)
#' - annots: List of annotations (annotation matrices)
#'
#' @examples
#' # load VDX dataset
#' data(vdxs)
#' # load NKI dataset
#' data(nkis)
#' # reduce datasets
#' ginter <- intersect(annot.vdxs[ ,"EntrezGene.ID"], annot.nkis[ ,"EntrezGene.ID"])
#' ginter <- ginter[!is.na(ginter)][1:30]
#' myx <- unique(c(match(ginter, annot.vdxs[ ,"EntrezGene.ID"]),
#'   sample(x=1:nrow(annot.vdxs), size=20)))
#' data2.vdxs <- data.vdxs[ ,myx]
#' annot2.vdxs <- annot.vdxs[myx, ]
#' myx <- unique(c(match(ginter, annot.nkis[ ,"EntrezGene.ID"]),
#'   sample(x=1:nrow(annot.nkis), size=20)))
#' data2.nkis <- data.nkis[ ,myx]
#' annot2.nkis <- annot.nkis[myx, ]
#' # mapping of datasets
#' datas <- list("VDX"=data2.vdxs,"NKI"=data2.nkis)
#' annots <- list("VDX"=annot2.vdxs, "NKI"=annot2.nkis)
#' datas.mapped <- map.datasets(datas=datas, annots=annots, do.mapping=TRUE)
#' str(datas.mapped, max.level=2)
#'
#' @md
#' @export
map.datasets <-
function(datas, annots, do.mapping=FALSE, mapping.coln="EntrezGene.ID", mapping, verbose=FALSE) {
	if((length(datas) != length(annots)) || !all(names(datas) == names(annots))) { stop("discordance between lists of datasets and annotations!") }
	## do the mapping (or not) and collect the set of unique features
	datas2 <- annots2 <- comid <- NULL
	for(k in 1:length(datas)) {
		if(verbose) { message(sprintf("%s", names(datas)[k])) }
		if(do.mapping) {
			gid <- as.character(annots[[k]][ ,mapping.coln])
			names(gid) <- dimnames(annots[[k]])[[1]]
			ugid <- unique(gid)
			ugid <- ugid[!is.na(ugid)]
			names(ugid) <- paste("geneid", ugid, sep=".")
			rr <- geneid.map(geneid1=gid, data1=datas[[k]], geneid2=ugid, verbose=FALSE)
			tt <- rr$data1
			## update gene ids since only missing values may be present for some of them
			ugid <- rr$geneid2
			dimnames(tt)[[2]] <- names(ugid)
			datas2 <- c(datas2, list(tt))
			tt <- annots[[k]][names(rr$geneid1), , drop=FALSE]
			dimnames(tt)[[1]] <- names(ugid)
			annots2 <- c(annots2, list(tt))
			comid <- unique(c(comid, names(ugid)))
			rm(rr)
			gc()
		} else {
			datas2 <- c(datas2, list(datas[[k]]))
			annots2 <- c(annots2, list(annots[[k]]))
			comid <- unique(c(comid, dimnames(datas[[k]])[[2]]))
		}
	}
	names(datas2) <- names(annots2) <- names(datas)
	#comid <- sort(comid)
	## put NA values for missing features
	for(k in 1:length(datas)) {
		tt <- matrix(NA, nrow=nrow(datas2[[k]]), ncol=length(comid), dimnames=list(dimnames(datas2[[k]])[[1]], comid))
		tt[dimnames(datas2[[k]])[[1]], dimnames(datas2[[k]])[[2]]] <- datas2[[k]]
		datas2[[k]] <- tt
		tt <- rbind(annots2[[k]], matrix(NA, nrow=length(comid) - nrow(annots2[[k]]), ncol=ncol(annots2[[k]]), dimnames=list(comid[!is.element(comid, dimnames(annots2[[k]])[[1]])], dimnames(annots2[[k]])[[2]])))
		tt <- tt[comid, , drop=FALSE]
		annots2[[k]] <- tt
	}
	return(list("datas"=datas2, "annots"=annots2))
}