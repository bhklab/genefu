\name{oncotypedx}
\alias{oncotypedx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the OncotypeDX signature as published by Paik et al. in 2004.
}
\description{
This function computes signature scores and risk classifications from gene expression values following the algorithm used for the OncotypeDX signature as published by Paik et al. 2004.
}
\usage{
oncotypedx(data, annot, do.mapping = FALSE, mapping, verbose = FALSE)
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
\details{
Note that for Affymetrix HGU datasets, the mapping is not necessary.
}
\value{
	\item{score }{Continuous signature scores}
	\item{risk }{Binary risk classification, \code{1} being high risk and \code{0} being low risk.}
	\item{mapping }{Mapping used if necessary.}
  	\item{probe }{If mapping is performed, this matrix contains the correspondence between the gene list (aka signature) and gene expression data.}
}
\references{
S. Paik, S. Shak, G. Tang, C. Kim, J. Bakker, M. Cronin, F. L. Baehner, M. G. Walker, D. Watson, T. Park, W. Hiller, E. R. Fisher, D. L. Wickerham, J. Bryant, and N. Wolmark (2004) "A Multigene Assay to Predict Recurrence of Tamoxifen-Treated, Node-Negative Breast Cancer", \emph{New England Journal of Medicine}, \bold{351}(27):2817--2826.
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
## load GENE70 signature
data(sig.oncotypedx)
## load NKI dataset
data(nki)
## compute relapse score
rs.nki <- oncotypedx(data=data.nki, annot=annot.nki, do.mapping=TRUE)
table(rs.nki$risk)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ prognosis }
%%\keyword{ breast cancer }
