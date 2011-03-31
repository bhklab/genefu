`map.datasets` <-
function(datas, annots, do.mapping=FALSE, mapping, verbose=FALSE) {
	if((length(datas) != length(annots)) || !all(names(datas) == names(annots))) { stop("discordance between lists of datasets and annotations!") }
	## do the mapping (or not) and collect the set of unique features
	datas2 <- annots2 <- comid <- NULL
	for(k in 1:length(datas)) {
		if(verbose) { message(sprintf("%s", names(datas)[k])) }
		if(do.mapping) {
			gid <- as.numeric(as.character(annots[[k]][ ,"EntrezGene.ID"]))
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