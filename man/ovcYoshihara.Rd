\name{ovcYoshihara}
\alias{ovcYoshihara}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the subtype scores and risk classifications for the prognostic signature published by Yoshihara et al.
}
\description{
This function computes subtype scores and risk classifications from gene expression values following the algorithm developed by Yoshihara et al, for prognosis in ovarian cancer.
}
\usage{
ovcYoshihara(data, annot, hgs, 
gmap = c("entrezgene", "ensembl_gene_id", "hgnc_symbol", "unigene", "refseq_mrna"), 
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
}
\references{
Yoshihara K, Tajima A, Yahata T, Kodama S, Fujiwara H, Suzuki M, Onishi Y, Hatae M, Sueyoshi K, Fujiwara H, Kudo, Yoshiki, Kotera K, Masuzaki H, Tashiro H, Katabuchi H, Inoue I, Tanaka K (2010) "Gene expression profile for predicting survival in advanced-stage serous ovarian cancer across two independent datasets", \emph{PloS one}, \bold{5}(3):e9615.
}
\author{
Benjamin Haibe-Kains
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[genefu]{sigOvcYoshihara}}
}
\examples{
## load the ovcYoshihara signature
data(sigOvcYoshihara)
## load NKI dataset
data(nkis)
colnames(annot.nkis)[is.element(colnames(annot.nkis), "EntrezGene.ID")] <- "entrezgene"
## compute relapse score
ovcYoshihara.nkis <- ovcYoshihara(data=data.nkis, 
annot=annot.nkis, gmap="entrezgene", do.mapping=TRUE)
table(ovcYoshihara.nkis$risk)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ prognosis }
%%\keyword{ breast cancer }
