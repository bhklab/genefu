#' @title Function to compute the Gene Expression progNostic Index Using Subtypes (GENIUS)
#'   as published by Haibe-Kains et al. 2010
#'
#' @description
#' This function computes the Gene Expression progNostic Index Using Subtypes (GENIUS)
#'   as published by Haibe-Kains et al. 2010. Subtype-specific risk scores are computed for
#'   each subtype signature separately and an overall risk score is computed by combining
#'   these scores with the posterior probability to belong to each of the breast cancer
#'   molecular subtypes.
#'
#' @usage
#' genius(data, annot, do.mapping = FALSE, mapping, do.scale = TRUE)
#'
#' @param data Matrix of gene expressions with samples in rows and probes in columns,
#'   dimnames being properly defined.
#' @param annot	Matrix of annotations with at least one column named "EntrezGene.ID",
#'   dimnames being properly defined.
#' @param do.mapping TRUE if the mapping through Entrez Gene ids must be performed (in case of ambiguities,
#'   the most variant probe is kept for each gene), FALSE otherwise.
#' @param mapping Matrix with columns "EntrezGene.ID" and "probe" used to force the
#'   mapping such that the probes are not selected based on their variance.
#' @param do.scale TRUE if the ESR1, ERBB2 and AURKA (module) scores must be rescaled
#'   (see rescale), FALSE otherwise.
#'
#' @return
#' A list with items:
#' - GENIUSM1: Risk score from the ER-/HER2- subtype signature in GENIUS model.
#' - GENIUSM2: Risk score from the HER2+ subtype signature in GENIUS model.
#' - GENIUSM3: Risk score from the ER+/HER2- subtype signature in GENIUS model.
#' - score: Overall risk prediction as computed by the GENIUS model.a.
#'
#' @references
#' Haibe-Kains B, Desmedt C, Rothe F, Sotiriou C and Bontempi G (2010) "A fuzzy gene
#' expression-based computational approach improves breast cancer prognostication", Genome Biology, 11(2):R18
#'
#' @seealso
#' [genefu::subtype.cluster.predict],[genefu::sig.score]
#'
#' @examples
#' # load NKI dataset
#' data(nkis)
#' # compute GENIUS risk scores based on GENIUS model fitted on VDX dataset
#' genius.nkis <- genius(data=data.nkis, annot=annot.nkis, do.mapping=TRUE)
#' str(genius.nkis)
#' # the performance of GENIUS overall risk score predictions are not optimal
#' # since only part of the NKI dataset was used
#'
#' @md
#' @export
genius <- function(data, annot, do.mapping=FALSE, mapping, do.scale=TRUE) {

	## predict breast cancer molecular subtypes
	sbt.id <- subtype.cluster.predict(sbt.model=scmod1.robust, data=data, annot=annot, do.mapping=do.mapping, mapping=mapping, do.prediction.strength=FALSE, do.BIC=FALSE, plot=FALSE, verbose=FALSE)

	usbt <- unique(sbt.id$subtype)
	usbt <- sort(usbt[!is.na(usbt)])
	pred.sbtclassif <- NULL
	for(ii in 1:length(usbt)) {
		myx <- sbt.id$subtype == usbt[ii] & !is.na(sbt.id$subtype)

		#compute the score from model
		score <- sig.score(x=sig.genius[[ii]][ , c("probe", "EntrezGene.ID",  "coefficient")], data=data, annot=annot, do.mapping=do.mapping,  mapping=mapping, verbose=FALSE)$score
		if(do.scale) {
			#the rescaling needs a large sample size!
			#necessary if we want to validate the classifier using a different dataset
			#the estimation of survival probabilities depends on the scale of the score
			score <-  (rescale(score, q=0.05, na.rm=TRUE) - 0.5) * 2
		}
		names(score) <- dimnames(data)[[1]]
		pred.sbtclassif <- c(pred.sbtclassif, list("score"=score))
	}
	names(pred.sbtclassif) <- names(sig.genius)
	#combine classifications
	cc <- NULL
	for(j in 1:length(pred.sbtclassif)) {
		cc <- cbind(cc, pred.sbtclassif[[j]])
	}
	ww <- sbt.id$subtype.proba
	combine.pred <- apply(ww * cc, 1, sum)

	pred.sbtclassif <- c(pred.sbtclassif, list(combine.pred))
	names(pred.sbtclassif)[length(pred.sbtclassif)] <- "score"
	return(pred.sbtclassif)
}