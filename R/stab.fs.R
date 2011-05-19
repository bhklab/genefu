`stab.fs` <-
function(fsets, N, method=c("kuncheva", "davis"), ...) {
	
	####################
	## internal functions
	####################
	
	kuncheva.stab <- function(fsets, N) {
		kk <- length(fsets)
		KI <- function(f1, f2, ss, NN) {
			#if(length(f1) != length(f2)) { stop("length of the two sets of selected features must be identical!") }
			#ss <- length(f1)
			if(ss == NN) { return(NA) }
			rr <- length(intersect(f1, f2))
			ki.est <- (rr - (ss^2 / NN)) / (ss - (ss^2 / NN))
			return(ki.est)
		}
		ss <- unique(unlist(lapply(fsets, length)))
		if(length(ss) > 1) { stop("length of sets of selected features must be identical!") }
		stab.res <- 0
		for(i in 1:(kk - 1)) {
			for(j in (i + 1):kk) {
				stab.res <- stab.res + KI(f1=fsets[[i]], f2=fsets[[j]], ss=ss, NN=N)
			}
		}
		return((2 * stab.res) / (kk * (kk - 1)))	
	}
	
	davis.stab <- function(fsets, N, penalty=1) {
		kk <- length(fsets)
		ss <- unique(unlist(lapply(fsets, length)))
		if(length(ss) > 1) { stop("length of sets of selected features must be identical!") }
		stab.res <- sum(sort(table(unlist(fsets)), decreasing=TRUE)[1:ss]) / (kk * ss)
		return(stab.res - penalty * (ss / N))	
	}
	
	####################
	
	method <- match.arg(method)
	if(!is.list(fsets)) { stop("fsets must be a list of sets of selected features!") }
	switch(method,
		"kuncheva"={
			stab <- kuncheva.stab(fsets=fsets, N=N)
		}, 
		"davis"={
			 stab <- davis.stab(fsets=fsets, N=N, ...)
		})
		return(stab)
}

## k <- 1000; fsets <- NULL; for(i in 1:k) { fsets <- c(fsets, list(sample(paste("feature", 1:10000, sep="."), 200))) }; names(fsets) <- paste("rand", 1:k, sep=".")
## kuncheva.stab(fsets=fsets, N=10000)
## davis.stab(fsets=fsets, N=10000, penalty=1)