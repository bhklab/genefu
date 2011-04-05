\name{compute.proto.cor.meta}
\alias{compute.proto.cor.meta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute correlations to prototypes in a meta-analytical framework
}
\description{
This function computes meta-estimate of correlation coefficients between a set of genes and a set of prototypes from a list of gene expression datasets.
}
\usage{
compute.proto.cor.meta(datas, proto, method = c("pearson", "spearman"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datas}{
List of datasets. Each dataset is a matrix of gene expressions with samples in rows and probes in columns, dimnames being properly defined. All the datasets must have the same probes.
}
  \item{proto}{
Names of prototypes (e.g. their EntrezGene ID).
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
  \item{cor }{Matrix of meta-estimate of correlation coefficients with probes in rows and prototypes in columns.}
  \item{cor.n }{Number of samples used to compute meta-estimate of correlation coefficients.}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Benjamin Haibe-Kains
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[genefu]{map.datasets}}
}
\examples{
## load VDX dataset
data(vdx)
## load NKI dataset
data(nki)
## reduce datasets
ginter <- intersect(annot.vdx[ ,"EntrezGene.ID"], annot.nki[ ,"EntrezGene.ID"])
ginter <- ginter[!is.na(ginter)][1:30]
myx <- unique(c(match(ginter, annot.vdx[ ,"EntrezGene.ID"]),
  sample(x=1:nrow(annot.vdx), size=20)))
data2.vdx <- data.vdx[ ,myx]
annot2.vdx <- annot.vdx[myx, ]
myx <- unique(c(match(ginter, annot.nki[ ,"EntrezGene.ID"]),
  sample(x=1:nrow(annot.nki), size=20)))
data2.nki <- data.nki[ ,myx]
annot2.nki <- annot.nki[myx, ]
## mapping of datasets
datas <- list("VDX"=data2.vdx,"NKI"=data2.nki)
annots <- list("VDX"=annot2.vdx, "NKI"=annot2.nki)
datas.mapped <- map.datasets(datas=datas, annots=annots, do.mapping=TRUE)
## define some prototypes
protos <- paste("geneid", ginter[1:3], sep=".")
## compute meta-estimate of correlation coefficients to the three prototype genes
probecor <- compute.proto.cor.meta(datas=datas.mapped$datas, proto=protos,
  method="pearson")
str(probecor)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ correlation }
