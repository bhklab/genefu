\name{bimod}
\alias{bimod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to identify bimodality for gene expression or signature score
}
\description{
This function fits a mixture of two Gaussians to identify bimodality. Useful to identify ER of HER2 status of breast tumors using ESR1 and ERBB2 expressions respectively.
}
\usage{
bimod(x, data, annot, do.mapping = FALSE, mapping, model = c("E", "V"), 
  do.scale = TRUE, verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Matrix containing the gene(s) in the gene list in rows and at least three columns: "probe", "EntrezGene.ID" and "coefficient" standing for the name of the probe, the NCBI Entrez Gene id and the coefficient giving the direction and the strength of the association of each gene in the gene list.
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
  \item{model}{
Model name used in \code{\link[mclust]{Mclust}}.
}
  \item{do.scale}{
\code{TRUE} if the gene expressions or signature scores must be rescaled (see \code{\link[genefu]{rescale}}), \code{FALSE} otherwise.
}
  \item{verbose}{
\code{TRUE} to print informative messages, \code{FALSE} otherwise.
}
  \item{\dots}{
Additional parameters to pass to \code{\link[genefu]{sig.score}}.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
  \item{status }{Status being 0 or 1.}
  \item{status1.proba }{Probability p to be of status 1, the probability to be of status 0 being 1-p.}
  \item{gaussians }{Matrix of parameters fitted in the mixture of two Gaussians. Matrix of NA values if EM algorithm did not converge.}
  \item{BIC }{Values (gene expressions or signature scores) used to identify bimodality.}
	\item{BI }{Bimodality Index (BI) as defined by Wang et al., 2009.}
  \item{x }{Values (gene expressions or signature scores) used to identify bimodality.}
}
\references{
Desmedt C, Haibe-Kains B, Wirapati P, Buyse M, Larsimont D, Bontempi G, Delorenzi M, Piccart M, and Sotiriou C (2008) "Biological processes associated with breast cancer clinical outcome depend on the molecular subtypes", \emph{Clinical Cancer Research}, \bold{14}(16):5158--5165.

Wirapati P, Sotiriou C, Kunkel S, Farmer P, Pradervand S, Haibe-Kains B, Desmedt C, Ignatiadis M, Sengstag T, Schutz F, Goldstein DR, Piccart MJ and Delorenzi M (2008) "Meta-analysis of Gene-Expression Profiles in Breast Cancer: Toward a Unified Understanding of Breast Cancer Sub-typing and Prognosis Signatures", \emph{Breast Cancer Research}, \bold{10}(4):R65.

Fraley C and Raftery E (2002) "Model-Based Clustering, Discriminant Analysis, and Density Estimation", \emph{Journal of American Statistical Asscoiation}, \bold{97}(458):611--631.

Wang J, Wen S, Symmans FW, Pusztai L and Coombes KR (2009) "The bimodality index: a criterion for discovering and ranking bimodal signatures from cancer gene expression profiling data", \emph{Cancer Informatics}, \bold{7}:199--216.
}
\author{
Benjamin Haibe-Kains
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[mclust]{Mclust}}
}
\examples{
## load NKI data
data(nki)
## load gene modules from Desmedt et al. 2008
data(mod1)
## retrieve esr1 affy probe and Entrez Gene id 
esr1 <- mod1$ESR1[1, ,drop=FALSE]
## computation of signature scores
esr1.bimod <- bimod(x=esr1, data=data.nki, annot=annot.nki, do.mapping=TRUE, 
  model="V", verbose=TRUE)
table("ER.IHC"=demo.nki[ ,"er"], "ER.GE"=esr1.bimod$status)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
