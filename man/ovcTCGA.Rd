\name{ovcTCGA}
\alias{ovcTCGA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the prediction scores and risk classifications for the ovarian cancer TCGA signature
}
\description{
This function computes signature scores and risk classifications from gene expression values following the algorithm developed by the TCGA consortium for ovarian cancer.
}
\usage{
ovcTCGA(data, annot, 
gmap = c("entrezgene", "ensembl_gene_id", "hgnc_symbol", "unigene"), 
do.mapping = FALSE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Matrix of gene expressions with samples in rows and probes in columns, dimnames being properly defined.
}
  \item{annot}{
Matrix of annotations with one column named as gmap, dimnames being properly defined.
}
  \item{gmap}{
character string containing the \code{biomaRt} attribute to use for mapping if \code{do.mapping=TRUE}
}
  \item{do.mapping}{
\code{TRUE} if the mapping through Entrez Gene ids must be performed (in case of ambiguities, the most variant probe is kept for each gene), \code{FALSE} otherwise.
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
Bell D, Berchuck A, Birrer M et al. (2011) "Integrated genomic analyses of ovarian carcinoma", \emph{Nature}, \bold{474}(7353):609--615
}
\author{
Benjamin Haibe-Kains
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[genefu]{sigOvcTCGA}}
}
\examples{
## load the ovcTCGA signature
data(sigOvcTCGA)
## load NKI dataset
data(nkis)
colnames(annot.nkis)[is.element(colnames(annot.nkis), "EntrezGene.ID")] <- "entrezgene"
## compute relapse score
ovcTCGA.nkis <- ovcTCGA(data=data.nkis, annot=annot.nkis, gmap="entrezgene", do.mapping=TRUE)
table(ovcTCGA.nkis$risk)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ prognosis }
%%\keyword{ breast cancer }
