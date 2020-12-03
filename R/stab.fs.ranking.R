#' @title Function to quantify stability of feature ranking
#'
#' @description
#' This function computes several indexes to quantify feature ranking 
#'   stability for several number of selected features. This is usually 
#'   estimated through perturbation of the original dataset by generating 
#'   multiple sets of selected features.
#'
#' @usage
#' stab.fs.ranking(fsets, sizes, N, method = c("kuncheva", "davis"), ...)
#' 
#' @param fsets	list or matrix of sets of selected features (in rows), 
#'   each ranking must have the same size.
#' @param sizes	Number of top-ranked features for which the stability 
#'   index must be computed.
#' @param N	total number of features on which feature selection is performed
#' @param method	stability index (see details section).
#' @param ...	additional parameters passed to stability index (penalty 
#'   that is a numeric for Davis' stability index, see details section).
#'
#' @details
#' Stability indices may use different parameters. In this version only the 
#'   Davis index requires an additional parameter that is penalty, a numeric 
#'   value used as penalty term.
#' Kuncheva index (kuncheva) lays in [-1, 1], An index of -1 means no 
#'   intersection between sets of selected features, +1 means that all the 
#'   same features are always selected and 0 is the expected stability of a 
#'   random feature selection.
#' Davis index (davis) lays in [0,1], With a pnalty term equal to 0, an index 
#'   of 0 means no intersection between sets of selected features and +1 means 
#'   that all the same features are always selected. A penalty of 1 is usually 
#'   used so that a feature selection performed with no or all features has a 
#'   Davis stability index equals to 0. None estimate of the expected Davis 
#'   stability index of a random feature selection was published.
#'
#' @return
#' A vector of numeric that are stability indices for each size of the sets 
#'   of selected features given the rankings.
#'
#' @references
#' Davis CA, Gerick F, Hintermair V, Friedel CC, Fundel K, Kuffner R, 
#'   Zimmer R (2006) "Reliable gene signatures for microarray classification: 
#'  assessment of stability and performance", Bioinformatics, 22(19):356-2363.
#' Kuncheva LI (2007) "A stability index for feature selection", AIAP'07: 
#'   Proceedings of the 25th conference on Proceedings of the 25th IASTED 
#'   International Multi-Conference, pages 390-395.
#'
#' @seealso
#' [genefu::stab.fs]
#' 
#' @examples
#' # 100 random selection of 50 features from a set of 10,000 features
#' fsets <- lapply(as.list(1:100), function(x, size=50, N=10000) {
#'   return(sample(1:N, size, replace=FALSE))} )
#' names(fsets) <- paste("fsel", 1:length(fsets), sep=".")
#' 
#' # Kuncheva index
#' stab.fs.ranking(fsets=fsets, sizes=c(1, 10, 20, 30, 40, 50),
#'   N=10000, method="kuncheva")
#' # close to 0 as expected for a random feature selection
#' 
#' # Davis index
#' stab.fs.ranking(fsets=fsets, sizes=c(1, 10, 20, 30, 40, 50),
#'   N=10000, method="davis", penalty=1)
#' 
#' @md
#' @export
stab.fs.ranking <-
function(fsets, sizes, N, method=c("kuncheva", "davis"), ...) {
	
	####################
	## internal functions
	####################
	
	kuncheva.stab.ranking <- function(fsets, N, x) {
		ss <- x
		fsets <- fsets[ , 1:ss, drop=FALSE]
		kk <- nrow(fsets)
		KI <- function(f1, f2, ss, NN) {
			#if(length(f1) != length(f2)) { stop("length of the two sets of selected features must be identical!") }
			#ss <- length(f1)
			if(ss == NN) { return(NA) }
			rr <- length(intersect(f1, f2))
			ki.est <- (rr - (ss^2 / NN)) / (ss - (ss^2 / NN))
			return(ki.est)
		}
		
		stab.res <- 0
		for(i in 1:(kk - 1)) {
			for(j in (i + 1):kk) {
				stab.res <- stab.res + KI(f1=fsets[i, ], f2=fsets[j, ], ss=ss, NN=N)
			}
		}
		return((2 * stab.res) / (kk * (kk - 1)))	
	}
	
	davis.stab.ranking <- function(fsets, N, x, penalty=1) {
		ss <- x
		fsets <- fsets[ , 1:ss, drop=FALSE]
		kk <- nrow(fsets)
		stab.res <- sum(sort(table(fsets), decreasing=TRUE)[1:ss]) / (kk * ss)
		return(stab.res - penalty * (ss / N))	
	}
	
	####################
	
	method <- match.arg(method)
	if(is.list(fsets)) { ## transform list into matrix
		Nn <- unique(unlist(lapply(fsets, length)))
		if(length(Nn) > 1) { stop("length of sets of selected features must be identical!") }
		nam <- names(fsets)
		fsets <- t(sapply(X=1:length(fsets), FUN=function(y, x) { return(y[[x]]) }, y=fsets))
		dimnames(fsets) <- list(nam, paste("rank", 1:Nn, sep="."))
	} else { Nn <- ncol(fsets) }
	if(missing(N)) { N <- Nn }
	if(missing(sizes)) { sizes <- 1:Nn }
	sizes <- sizes[sizes <= Nn]
	
	switch(method,
		"kuncheva"={
			stab <- unlist(sapply(X=sizes, FUN=kuncheva.stab.ranking, fsets=fsets, N=N))
		}, 
		"davis"={
			 stab <- unlist(sapply(X=sizes, FUN=davis.stab.ranking, fsets=fsets, N=N, ...))
		})
		names(stab) <- paste("size", sizes, sep=".")
		return(stab)
}

## k <- 100; fsets <- NULL; for(i in 1:k) { fsets <- c(fsets, list(sample(paste("feature", 1:1000, sep=".")))) }; names(fsets) <- paste("rand", 1:k, sep=".")
## stab.fs.ranking(fsets=fsets, sizes=1:10, method="kuncheva", penalty=1)