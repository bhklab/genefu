`gene76` <-
function(data, er) {
	
	A <- 313.5
	B <- 280
	
	score <- NULL
	for(i in 1:nrow(data)) {
		if(is.na(er[i])) { score <- c(score, NA) }
		else {
			if(er[i] == 1) {
				score <- c(score, A + sum(data[i,dimnames(sig.gene76)[[1]][sig.gene76[ ,"er"] == 1]] * sig.gene76[sig.gene76[ ,"er"] == 1,"std.cox.coefficient"]))
			}
			else {
				score <- c(score, B + sum(data[i,dimnames(sig.gene76)[[1]][sig.gene76[ ,"er"] == 0]] * sig.gene76[sig.gene76[ ,"er"] == 0,"std.cox.coefficient"]))
			}
		}
	}
	names(score) <- dimnames(data)[[1]]
	risk <- ifelse(score >= 0, 1, 0)	
	
	return(list("score"=score, "risk"=risk))
}
