\name{ggi}
\alias{ggi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the raw and scaled Gene expression Grade Index (GGI)
}
\description{
This function computes signature scores and risk classifications from gene expression values following the algorithm used for the Gene expression Grade Index (GGI).
}
\usage{
ggi(data, annot, do.mapping = FALSE, mapping, hg, verbose = FALSE)
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
  \item{hg}{
Vector containing the histological grade (HG) status of breast cancer patients in the dataset.
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
Sotiriou C, Wirapati P, Loi S, Harris A, Bergh J, Smeds J, Farmer P, Praz V, Haibe-Kains B, Lallemand F, Buyse M, Piccart MJ and Delorenzi M (2006) "Gene expression profiling in breast cancer: Understanding the molecular basis of histologic grade to improve prognosis", \emph{Journal of National Cancer Institute}, \bold{98}:262--272
}
\author{
Benjamin Haibe-Kains
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[genefu]{gene76}}
}
\examples{
## load GGI signature
data(sig.ggi)
## load NKI dataset
data(nki)
## compute relapse score
ggi.nki <- ggi(data=data.nki, annot=annot.nki, do.mapping=TRUE,
  hg=demo.nki[ ,"grade"])
table(ggi.nki$risk)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ prognosis }
%%\keyword{ breast cancer }
