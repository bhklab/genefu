\name{rorS}
\alias{rorS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the rorS signature as published by Parker et al 2009
}
\description{
This function computes signature scores and risk classifications from gene expression values following the algorithm used for the rorS signature as published by Parker et al 2009.
}
\usage{
rorS(data, annot, do.mapping = FALSE, mapping, verbose = FALSE)
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
\code{TRUE} if the mapping through Entrez Gene ids must be performed (in case of ambiguities, the most variant probe is kept for each gene), \code{FALSE} otherwise. Note that for Affymetrix HGU datasets, the mapping is not necessary.
}
  \item{mapping}{
Matrix with columns "EntrezGene.ID" and "probe" used to force the mapping such that the probes are not selected based on their variance.
}
  \item{verbose}{
\code{TRUE} to print informative messages, \code{FALSE} otherwise.
}
}
\value{
	\item{score }{Continuous signature scores}
	\item{risk }{Binary risk classification, \code{1} being high risk and \code{0} being low risk.}
	\item{mapping }{Mapping used if necessary.}
  	\item{probe }{If mapping is performed, this matrix contains the correspondence between the gene list (aka signature) and gene expression data.}
}
\references{
Parker, Joel S. and Mullins, Michael and Cheang, Maggie C.U. and Leung, Samuel and Voduc, David and Vickery, Tammi and Davies, Sherri and Fauron, Christiane and He, Xiaping and Hu, Zhiyuan and Quackenbush, John F. and Stijleman, Inge J. and Palazzo, Juan and Marron, J.S. and Nobel, Andrew B. and Mardis, Elaine and Nielsen, Torsten O. and Ellis, Matthew J. and Perou, Charles M. and Bernard, Philip S. (2009) "Supervised Risk Predictor of Breast Cancer Based on Intrinsic Subtypes", \emph{Journal of Clinical Oncology}, \bold{27}(8):1160--1167
}
\author{
Benjamin Haibe-Kains
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%%}

\examples{
## load NKI dataset
data(vdxs)
## compute relapse score
rs.vdxs <- rorS(data=data.vdxs, annot=annot.vdxs, do.mapping=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ prognosis }
%%\keyword{ breast cancer }
