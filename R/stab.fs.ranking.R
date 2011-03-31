`stab.fs.ranking` <-
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