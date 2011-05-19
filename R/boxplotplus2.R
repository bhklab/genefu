`boxplotplus2` <-
function(x, .jit = 0.25, .las = 1, .ylim, box.col = "lightgrey", pt.col = "blue", pt.cex=0.5, pt.pch=16, med.line = FALSE, med.col = "goldenrod", ...) {

    isMAT <- is.matrix(x)
    y <- x
    if (isMAT) {
		y <- data.frame(t(x))
		if(missing(.ylim)) { myrange <- range(y, na.rm=TRUE) } else { myrange <- .ylim }
	} else { if(missing(.ylim)) { myrange <- range(unlist(x), na.rm=TRUE) } else { myrange <- .ylim } }

    bp <- boxplot(y, las = .las, cex.axis = 0.85, border="grey", col=box.col, boxwex=0.5, ylim=myrange, range=0, ...)
    if (isMAT) {
        xp <- rep(1:nrow(x), times=ncol(x))
        yp <- as.vector(x)
    } else {
        reps <- sapply(x, FUN=function(x) length(x) )
        xp <- rep(1:length(y), times=reps)
        yp <- unlist(y)
    }
    points(jitter(xp, .jit), yp, cex=pt.cex, pch=pt.pch, col=pt.col)
    if (med.line) points(1:length(bp$n), bp$stats[3, ], type="b", col="goldenrod", lwd=3, pch=19)
    n <- bp$n
    names(n) <- bp$names
    n
}