\name{ovcAngiogenic}
\alias{ovcAngiogenic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the subtype scores and risk classifications for the angiogenic molecular subtype in ovarian cancer
}
\description{
This function computes subtype scores and risk classifications from gene expression values following the algorithm developed by Bentink, Haibe-Kains et al. to identify the angiogenic molecular subtype in ovarian cancer.
}
\usage{
ovcAngiogenic(data, annot, hgs, gmap = c("entrezgene", "ensembl_gene_id", "hgnc_symbol", "unigene"), do.mapping = FALSE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Matrix of gene expressions with samples in rows and probes in columns, dimnames being properly defined.
}
  \item{annot}{
Matrix of annotations with at least one column named "EntrezGene.ID", dimnames being properly defined.
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
	\item{subtype }{data frame reporting the subtype score, maximum likelihood classification and corresponding subtype probabilities}
}
\references{
Bentink S, Haibe-Kains B, Risch T, Fan J-B, Hirsch MS, Holton K, Rubio R, April C, Chen J, Wickham-Garcia E, Liu J, Culhane AC, Drapkin R, Quackenbush JF, Matulonis UA (2012) "Angiogenic mRNA and microRNA Gene Expression Signature Predicts a Novel Subtype of Serous Ovarian Cancer", \emph{PloS one}, \bold{7}(2):e30269

}
\author{
Benjamin Haibe-Kains
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[genefu]{sigOvcAngiogenic}}
}
\examples{
## load GGI signature
data(sigOvcAngiogenic)
## load NKI dataset
data(nkis)
colnames(annot.nkis)[colnames(annot.nkis) %in% "EntrezGene.ID"] <- "entrezgene"
## compute relapse score
ovcAngiogenic.nkis <- ovcAngiogenic(data=data.nkis, annot=annot.nkis, gmap="entrezgene", do.mapping=TRUE)
table(ovcAngiogenic.nkis$risk)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ prognosis }
%%\keyword{ breast cancer }
