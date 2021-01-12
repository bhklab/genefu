#' @title Box plot of group of values with corresponding jittered points
#'
#' @description
#' This function allows for display a boxplot with jittered points.
#'
#' @usage
#' boxplotplus2(x, .jit = 0.25, .las = 1, .ylim, box.col = "lightgrey",
#'  pt.col = "blue", pt.cex = 0.5, pt.pch = 16, med.line = FALSE,
#'  med.col = "goldenrod", ...)
#'
#' @param x could be a list of group values or a matrix (each group is a row).
#' @param .jit Amount of jittering noise.
#' @param .las Numeric in 0,1,2,3; the style of axis labels.
#' @param .ylim Range for y axis.
#' @param box.col Color for boxes.
#' @param pt.col Color for groups (jittered points).
#' @param pt.cex A numerical value giving the amount by which plotting jittered points
#'   should be magnified relative to the default.
#' @param pt.pch Either an integer specifying a symbol or a single character to be used
#'   as the default in plotting jittered points. See points for possible values and
#'   their interpretation.
#' @param med.line TRUE if a line should link the median of each group, FALSE otherwise.
#' @param med.col Color of med.line.
#' @param ... Additional parameters for boxplot function.
#'
#' @return
#' Number of samples in each group.
#'
#' @note
#' 2.21.2006 - Christos Hatzis, Nuvera Biosciences
#'
#' @seealso
#' [graphics::boxplot], [base::jitter]
#'
#' @examples
#' dd <- list("G1"=runif(20), "G2"=rexp(30) * -1.1, "G3"=rnorm(15) * 1.3)
#' boxplotplus2(x=dd, .las=3, .jit=0.75, .ylim=c(-3,3), pt.cex=0.75,
#'   pt.col=c(rep("darkred", 20), rep("darkgreen", 30), rep("darkblue", 15)),
#'   pt.pch=c(0, 9, 17))
#'
#' @md
#' @importFrom survcomp fisherz
#' @importFrom graphics points
#' @export
boxplotplus2 <- function(x, .jit = 0.25, .las = 1, .ylim,
    box.col="lightgrey", pt.col="blue", pt.cex=0.5, pt.pch=16,
    med.line = FALSE, med.col = "goldenrod", ...)
{

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