\name{ovcCrijns}
\alias{ovcCrijns}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the subtype scores and risk classifications for the prognostic signature published by Crinjs et al.
}
\description{
This function computes subtype scores and risk classifications from gene expression values using teh weights published by Crijns et al.
}
\usage{
ovcCrijns(data, annot, hgs, gmap = c("entrezgene", "ensembl_gene_id", "hgnc_symbol", "unigene"), do.mapping = FALSE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Matrix of gene expressions with samples in rows and probes in columns, dimnames being properly defined.
}
  \item{annot}{
Matrix of annotations with one column named as gmap, dimnames being properly defined.
}
  \item{hgs}{
vector of booleans with TRUE represents the ovarian cancer patients who have a high grade, late stage, serous tumor, FALSE otherwise. This is particularly important for properly rescaling the data. If \code{hgs} is missing, all the patients will be used to rescale the subtype score.
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
\details{
Note that the original algorithm has not been implemented as it necessitates refitting of the model weights in each new dataset. However the current implementation should give similar results.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
	\item{score }{Continuous signature scores}
	\item{risk }{Binary risk classification, \code{1} being high risk and \code{0} being low risk.}
	\item{mapping }{Mapping used if necessary.}
	\item{probe }{If mapping is performed, this matrix contains the correspondence between the gene list (aka signature) and gene expression data.}
}
\references{
Crijns APG, Fehrmann RSN, de Jong S, Gerbens F, Meersma G J, Klip HG, Hollema H, Hofstra RMW, te Meerman GJ, de Vries EGE, van der Zee AGJ (2009) "Survival-Related Profile, Pathways, and Transcription Factors in Ovarian Cancer"
\emph{PLoS Medicine}, \bold{6}(2):e1000024.
}
\author{
Benjamin Haibe-Kains
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[genefu]{sigOvcCrijns}}
}
\examples{
## load the ovsCrijns signature
data(sigOvcCrijns)
## load NKI dataset
data(nkis)
colnames(annot.nkis)[is.element(colnames(annot.nkis), "EntrezGene.ID")] <- "entrezgene"
## compute relapse score
ovcCrijns.nkis <- ovcCrijns(data=data.nkis, annot=annot.nkis, gmap="entrezgene", do.mapping=TRUE)
table(ovcCrijns.nkis$risk)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ prognosis }
%%\keyword{ breast cancer }
