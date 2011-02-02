`write.m.file` <-
function(obj, file, ...) {
	lcn <- dimnames(obj[[1]])[[2]]
	c1 <- c2 <- NULL
	for (i in 1:length(obj)) {
		ct <- names(obj)[i]
		tt <- as.matrix(obj[[i]])
		if(ncol(tt) == 1) { tt <- t(tt) }
		c1 <- c(c1, ct, rep("", nrow(tt)))
		c2 <- rbind(c2, obj[[i]], rep("", ncol(tt)))
	}
	dimnames(c2)[[1]] <- 1:nrow(c2)
	res <- cbind(c1, c2)[-length(c1), ,drop=FALSE]
	dimnames(res)[[2]] <- c("gene.list", lcn)
	dimnames(res)[[1]] <- 1:nrow(res)
	write.table(res, file=file, row.names=FALSE, sep=",", ...)
	invisible(res)
}