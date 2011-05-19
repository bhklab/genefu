\name{subtype.cluster.predict}
\alias{subtype.cluster.predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to identify breast cancer molecular subtypes using the Subtype Clustering Model
}
\description{
This function identifies the breast cancer molecular subtypes using a Subtype Clustering Model fitted by \code{\link[genefu]{subtype.cluster}}.
}
\usage{
subtype.cluster.predict(sbt.model, data, annot, do.mapping = FALSE,
  mapping, do.scale = TRUE, do.prediction.strength = FALSE,
  do.BIC = FALSE, plot = FALSE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sbt.model}{
Subtype Clustering Model as returned by \code{\link[genefu]{subtype.cluster}}.
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
  \item{do.prediction.strength}{
\code{TRUE} if the prediction strength must be computed (Tibshirani and Walther 2005), \code{FALSE} otherwise.
}
  \item{do.BIC}{
\code{TRUE} if the Bayesian Information Criterion must be computed for number of clusters ranging from 1 to 10, \code{FALSE} otherwise.
}
  \item{plot}{
\code{TRUE} if the patients and their corresponding subtypes must be plotted, \code{FALSE} otherwise.
}
  \item{verbose}{
\code{TRUE} to print informative messages, \code{FALSE} otherwise.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
  \item{subtype }{Subtypes identified by the Subtype Clustering Model. Subtypes can be either "ER-/HER2-", "HER2+" or "ER+/HER2-".}
  \item{subtype.proba }{Probabilities  to belong to each subtype estimated by the Subtype Clustering Model.}
  \item{prediction.strength }{Prediction strength for \code{subtypes}.}
  \item{BIC }{Bayesian Information Criterion for the Subtype Clustering Model with number of clusters ranging from 1 to 10.}
  \item{subtype2 }{Subtypes identified by the Subtype Clustering Model using AURKA to discriminate low and high proliferative tumors. Subtypes can be either "ER-/HER2-", "HER2+", "ER+/HER2- High Prolif" or "ER+/HER2- Low Prolif".}
  \item{prediction.strength2 }{Prediction strength for \code{subtypes2}.}
  \item{module.scores }{Matrix containing ESR1, ERBB2 and AURKA module scores.}
  \item{mapping }{Mapping if necessary (list of matrices with 3 columns: probe, EntrezGene.ID and new.probe).}
}
\references{
Desmedt C, Haibe-Kains B, Wirapati P, Buyse M, Larsimont D, Bontempi G, Delorenzi M, Piccart M, and Sotiriou C (2008) "Biological processes associated with breast cancer clinical outcome depend on the molecular subtypes", \emph{Clinical Cancer Research}, \bold{14}(16):5158--5165.

Wirapati P, Sotiriou C, Kunkel S, Farmer P, Pradervand S, Haibe-Kains B, Desmedt C, Ignatiadis M, Sengstag T, Schutz F, Goldstein DR, Piccart MJ and Delorenzi M (2008) "Meta-analysis of Gene-Expression Profiles in Breast Cancer: Toward a Unified Understanding of Breast Cancer Sub-typing and Prognosis Signatures", \emph{Breast Cancer Research}, \bold{10}(4):R65.

Tibshirani R and Walther G (2005) "Cluster Validation by Prediction Strength", \emph{Journal of Computational and Graphical Statistics}, \bold{14}(3):511--528
}
\author{
Benjamin Haibe-Kains
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[genefu]{subtype.cluster}}, \code{\link[genefu]{scmod1}}, \code{\link[genefu]{scmod2}}
}
\examples{
## load Subtype Clustering Model fitted on VDX
data(scmod1)
## load NKI data
data(nki)
## Subtype Clustering Model fitted on VDX and applied on NKI
scmod1.nki <- subtype.cluster.predict(sbt.model=scmod1, data=data.nki,
  annot=annot.nki, do.mapping=TRUE, do.scale=TRUE,
  do.prediction.strength=FALSE, do.BIC=FALSE, plot=TRUE, verbose=TRUE)
table(scmod1.nki$subtype)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ clustering }
