`pik3cags` <-
function(data, annot, do.mapping=FALSE, mapping, verbose=FALSE) {

	pik3cags.gl <- sig.pik3cags[ ,c("probe", "EntrezGene.ID", "coefficient")]
	res <- sig.score(x=pik3cags.gl, data=data, annot=annot, do.mapping=do.mapping, mapping=mapping, signed=TRUE, verbose=verbose)$score
	
	return (res)
}
