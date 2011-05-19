\name{pik3cags}
\alias{pik3cags}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the PIK3CA gene signature (PIK3CA-GS)
}
\description{
This function computes signature scores from gene expression values following the algorithm used for the PIK3CA gene signature (PIK3CA-GS).
}
\usage{
pik3cags(data, annot, do.mapping = FALSE, mapping, verbose = FALSE)
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
Vector of signature scores for PIK3CA-GS
}
\references{
Loi S, Haibe-Kains B, Majjaj S, Lallemand F, Durbecq V, Larsimont D, Gonzalez-Angulo AM, Pusztai L, Symmans FW, Bardelli A, Ellis P, Tutt AN, Gillett CE, Hennessy BT., Mills GB, Phillips WA, Piccart MJ, Speed TP, McArthur GA, Sotiriou C (2010) "PIK3CA mutations associated with gene signature of low mTORC1 signaling and better outcomes in estrogen receptor-positive breast cancer", \emph{Proceedings of the National Academy of Sciences}, \bold{107}(22):10208--10213
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
data(sig.pik3cags)
## load NKI dataset
data(nki)
## compute relapse score
pik3cags.nki <- pik3cags(data=data.nki, annot=annot.nki, do.mapping=TRUE)
head(pik3cags.nki)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ prognosis }
%%\keyword{ breast cancer }
