\name{tamr13}
\alias{tamr13}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the risk scores of the tamoxifen resistance signature (TAMR13)
}
\description{
This function computes signature scores  from gene expression values following the algorithm used for the Tamoxifen Resistance signature (TAMR13).
}
\usage{
tamr13(data, annot, do.mapping = FALSE, mapping, verbose = FALSE)
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
  \item{risk }{Binary risk classification, \code{1} being high risk and \code{0} being low risk (not implemented, the function will return NA values).}
}
\references{
Loi S, Haibe-Kains B, Desmedt C, Wirapati P, Lallemand F, Tutt AM, Gillet C, Ellis P, Ryder K, Reid JF, Daidone MG, Pierotti MA, Berns EMJJ, Jansen MPHM, Foekens JA, Delorenzi M, Bontempi G, Piccart MJ and Sotiriou C (2008) "Predicting prognosis using molecular profiling in estrogen receptor-positive breast cancer treated with tamoxifen", \emph{BMC Genomics}, \bold{9}(1):239
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
## load TAMR13 signature
data(sig.tamr13)
## load VDX dataset
data(vdx)
## compute relapse score
tamr13.vdx <- tamr13(data=data.vdx, annot=annot.vdx, do.mapping=FALSE)
summary(tamr13.vdx$score)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ prognosis }
%%\keyword{ breast cancer }
