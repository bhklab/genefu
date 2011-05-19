\name{subtype.cluster}
\alias{subtype.cluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to fit the Subtype Clustering Model
}
\description{
This function fits the Subtype Clustering Model as published in Desmedt et al. 2008 and Wiarapati et al. 2008. This model is actually a mixture of three Gaussians with equal shape, volume and variance (see \code{EEI} model in \code{\link[mclust]{Mclust}}). This model is adapted to breast cancer and uses ESR1, ERBB2 and AURKA dimensions to identify the molecular subtypes, i.e. ER-/HER2-, HER2+ and ER+/HER2- (Low and High Prolif).
}
\usage{
subtype.cluster(module.ESR1, module.ERBB2, module.AURKA, data, annot,
  do.mapping = FALSE, mapping, do.scale = TRUE, rescale.q = 0.05,
  model.name = "EEI", do.BIC = FALSE, plot = FALSE, filen, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{module.ESR1}{
Matrix containing the ESR1-related gene(s) in rows and at least three columns: "probe", "EntrezGene.ID" and "coefficient" standing for the name of the probe, the NCBI Entrez Gene id and the coefficient giving the direction and the strength of the association of each gene in the gene list.
}
  \item{module.ERBB2}{
Idem for ERBB2.
}
  \item{module.AURKA}{
Idem for AURKA.
}
  \item{data}{
Matrix of gene expressions with samples in rows and probes in columns, dimnames being properly defined.
}
  \item{annot}{
Matrix of annotations with at least one column named "EntrezGene.ID", dimnames being properly defined.
}
  \item{do.mapping}{
\code{TRUE} if the mapping through Entrez Gene ids must be performed (in case of ambiguities, the most variant probe is kept for each gene), \code{FALSE} otherwise.
}
  \item{mapping}{
Matrix with columns "EntrezGene.ID" and "probe" used to force the mapping such that the probes are not selected based on  their variance.
}
  \item{do.scale}{
\code{TRUE} if the ESR1, ERBB2 and AURKA (module) scores must be rescaled (see \code{\link[genefu]{rescale}}), \code{FALSE} otherwise.
}
  \item{rescale.q}{
	Proportion of expected outliers for rescaling the gene expressions.
}
  \item{do.BIC}{
\code{TRUE} if the Bayesian Information Criterion must be computed for number of clusters ranging from 1 to 10, \code{FALSE} otherwise.
}
  \item{model.name}{
Name of the model used to fit the mixture of Gaussians with the \link[mclust]{Mclust} from the \code{mclust} package; default is "EEI" for fitting a mixture of Gaussians with diagonal variance, equal volume, equal shape and identical orientation.
}
  \item{plot}{
\code{TRUE} if the patients and their corresponding subtypes must be plotted, \code{FALSE} otherwise.
}
  \item{filen}{
Name of the csv file where the subtype clustering model must be stored.
}
  \item{verbose}{
\code{TRUE} to print informative messages, \code{FALSE} otherwise.
}
}
%%\details{
%%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{model }{Subtype Clustering Model (mixture of three Gaussians), like \code{\link[genefu]{scmgene}}, \code{\link[genefu]{scmod1}} and \code{\link[genefu]{scmod2}} when this function is used on VDX dataset (Wang et al. 2005 and Minn et al. 2007) with the gene modules published in the two references cited below.}
  \item{BIC }{Bayesian Information Criterion for the Subtype Clustering Model with number of clusters ranging from 1 to 10.}
  \item{subtype }{Subtypes identified by the Subtype Clustering Model. Subtypes can be either "ER-/HER2-", "HER2+" or "ER+/HER2-".}
  \item{subtype.proba }{Probabilities  to belong to each subtype estimated by the Subtype Clustering Model.}
  \item{subtype2 }{Subtypes identified by the Subtype Clustering Model using AURKA to discriminate low and high proliferative tumors. Subtypes can be either "ER-/HER2-", "HER2+", "ER+/HER2- High Prolif" or "ER+/HER2- Low Prolif".}
  \item{module.scores }{Matrix containing ESR1, ERBB2 and AURKA module scores.}
}
\references{
Desmedt C, Haibe-Kains B, Wirapati P, Buyse M, Larsimont D, Bontempi G, Delorenzi M, Piccart M, and Sotiriou C (2008) "Biological processes associated with breast cancer clinical outcome depend on the molecular subtypes", \emph{Clinical Cancer Research}, \bold{14}(16):5158--5165.

Wirapati P, Sotiriou C, Kunkel S, Farmer P, Pradervand S, Haibe-Kains B, Desmedt C, Ignatiadis M, Sengstag T, Schutz F, Goldstein DR, Piccart MJ and Delorenzi M (2008) "Meta-analysis of Gene-Expression Profiles in Breast Cancer: Toward a Unified Understanding of Breast Cancer Sub-typing and Prognosis Signatures", \emph{Breast Cancer Research}, \bold{10}(4):R65.
}
\author{
Benjamin Haibe-Kains
}
%%\note{
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[genefu]{subtype.cluster.predict}}, \code{\link[genefu]{intrinsic.cluster}},  \code{\link[genefu]{intrinsic.cluster.predict}}, \code{\link[genefu]{scmod1}}, \code{\link[genefu]{scmod2}}
}
\examples{
## example without gene mapping
## load VDX data
data(vdx)
## load gene modules
data(mod1)
## fit a Subtype Clustering Model
scmod1.vdx <- subtype.cluster(module.ESR1=mod1$ESR1, module.ERBB2=mod1$ERBB2,
  module.AURKA=mod1$AURKA, data=data.vdx, annot=annot.vdx, do.mapping=FALSE,
  do.scale=TRUE, plot=TRUE, verbose=TRUE)
str(scmod1.vdx, max.level=1)
table(scmod1.vdx$subtype2)
## example with gene mapping
## load NKI data
data(nki)
## load gene modules
data(mod1)
## fit a Subtype Clustering Model
scmod1.nki <- subtype.cluster(module.ESR1=mod1$ESR1, module.ERBB2=mod1$ERBB2,
  module.AURKA=mod1$AURKA, data=data.nki, annot=annot.nki, do.mapping=TRUE,
  do.scale=TRUE, plot=TRUE, verbose=TRUE)
str(scmod1.nki, max.level=1)
table(scmod1.nki$subtype2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ clustering }
