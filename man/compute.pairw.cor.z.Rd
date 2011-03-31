\name{compute.pairw.cor.z}
\alias{compute.pairw.cor.z}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the Z transformation of the pairwise correlations for a list of datasets
}
\description{
This function computes the Z transformation of the meta-estimate of pairwise correlation coefficients for a set of genes  from a list of gene expression datasets.
}
\usage{
compute.pairw.cor.z(datas, method = c("pearson"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{datas}{
List of datasets. Each dataset is a matrix of gene expressions with samples in rows and probes in columns, dimnames being properly defined. All the datasets must have the same probes.
}
  \item{method}{
Estimator for correlation coefficient, can be either \code{pearson} or \code{spearman}.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{z }{Z transformation of the meta-estimate of correlation coefficients.}
  \item{se }{Standard error of the Z transformation of the meta-estimate of correlation coefficients.}
  \item{nn }{Number of samples used to compute the meta-estimate of correlation coefficients.}
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Benjamin Haibe-Kains
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[genefu]{map.datasets}}, \code{\link[genefu]{compute.pairw.cor.meta}}, \code{\link[genefu]{compute.proto.cor.meta}}
}
%%\examples{
%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ correlation }
