\name{gene70}
\alias{gene70}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the 70 genes prognosis profile (GENE70) as published by van't Veer et al. 2002
}
\description{
This function computes signature scores and risk classifications from gene expression values following the algorithm used for the 70 genes prognosis profile (GENE70) as published by van't Veer et al. 2002.
}
\usage{
gene70(data, annot, do.mapping = FALSE, mapping,
  std = c("none", "scale", "robust"), verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
Matrix with columns "EntrezGene.ID" and "probe" used to force the mapping such that the probes are not selected based on their variance.
}
  \item{std}{
Standardization of gene expressions: \code{scale} for traditional standardization based on mean and standard deviation, \code{robust} for standardization based on the 0.025 and 0.975 quantiles, \code{none} to keep gene expressions unchanged.
}
  \item{verbose}{
\code{TRUE} to print informative messages, \code{FALSE} otherwise.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
	\item{score }{Continuous signature scores}
	\item{risk }{Binary risk classification, \code{1} being high risk and \code{0} being low risk.}
	\item{mapping }{Mapping used if necessary.}
  	\item{probe }{If mapping is performed, this matrix contains the correspondence between the gene list (aka signature) and gene expression data.}
}
\references{
L. J. van't Veer and H. Dai and M. J. van de Vijver and Y. D. He and A. A. Hart and M. Mao and H. L. Peterse and K. van der Kooy and M. J. Marton and A. T. Witteveen and G. J. Schreiber and R. M. Kerkhiven and C. Roberts and P. S. Linsley and R. Bernards and S. H. Friend (2002) "Gene Expression Profiling Predicts Clinical Outcome of Breast Cancer", \emph{Nature}, \bold{415}:530--536.
}
\author{
Benjamin Haibe-Kains
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[genefu]{nki}}
}
\examples{
## load GENE70 signature
data(sig.gene70)
## load NKI dataset
data(nki)
## compute relapse score
rs.nki <- gene70(data=data.nki)
table(rs.nki$risk)
## note that the discrepancies compared to the original publication
## are closed to the official cutoff, raising doubts on its exact value.
## computation of the signature scores on a different microarray platform
## load VDX dataset
data(vdx)
## compute relapse score
rs.vdx <- gene70(data=data.vdx, annot=annot.vdx, do.mapping=TRUE)
table(rs.vdx$risk)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ prognosis }
%%\keyword{ breast cancer }
