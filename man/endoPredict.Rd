\name{endoPredict}
\alias{endoPredict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the endoPredict signature as published by Filipits et al 2011
}
\description{
This function computes signature scores and risk classifications from gene expression values following the algorithm used for the endoPredict signature as published by Filipits et al 2011.
}
\usage{
endoPredict(data, annot, do.mapping = FALSE, mapping, verbose = FALSE)
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
The function works best if data have been noralized with MAS5.
Note that for Affymetrix HGU datasets, the mapping is not necessary.
}
\value{
	\item{score }{Continuous signature scores}
	\item{risk }{Binary risk classification, \code{1} being high risk and \code{0} being low risk.}
	\item{mapping }{Mapping used if necessary.}
  	\item{probe }{If mapping is performed, this matrix contains the correspondence between the gene list (aka signature) and gene expression data.}
}
\references{
ilipits, M., Rudas, M., Jakesz, R., Dubsky, P., Fitzal, F., Singer, C. F., et al. (2011). "A new molecular predictor of distant recurrence in ER-positive, HER2-negative breast cancer adds independent information to conventional clinical risk factors." \emph{Clinical Cancer Research}, \bold{17}(18):6012-–6020.
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
data(sig.endoPredict)
## load NKI dataset
data(vdxs)
## compute relapse score
rs.vdxs <- endoPredict(data=data.vdxs, annot=annot.vdxs, do.mapping=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ prognosis }
%%\keyword{ breast cancer }
