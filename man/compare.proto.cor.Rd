\name{compare.proto.cor}
\alias{compare.proto.cor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to statistically compare correlation to prototypes
}
\description{
This function performs a statistical comparison of the correlation coefficients as computed between each probe and prototype.
}
\usage{
compare.proto.cor(gene.cor, proto.cor, nn,
  p.adjust.m = c("none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene.cor}{
Correlation coefficients between the probes and each of the prototypes.
}
  \item{proto.cor}{
Pairwise correlation coefficients of the prototypes.
}
  \item{nn}{
Number of samples used to compute the correlation coefficients between the probes and each of the prototypes.
}
  \item{p.adjust.m}{
Correction method as defined in \code{\link[stats]{p.adjust}}.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
Data frame with probes in rows and with three columns: "proto" is the prototype to which the probe is the most correlated, "cor" is the actual correlation, and "signif" is the (corrected) p-value for the superiority of the correlation to this prototype compared to the second highest correlation.
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Benjamin Haibe-Kains
}
%%\note{
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[genefu]{compute.proto.cor.meta}}, \code{\link[genefu]{compute.pairw.cor.meta}}
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
## compute meta-estimate of pairwise correlation coefficients between prototypes
datas.proto <- lapply(X=datas.mapped$datas, FUN=function(x, p) {
  return(x[ ,p,drop=FALSE]) }, p=protos)
protocor <- compute.pairw.cor.meta(datas=datas.proto, method="pearson")
## compare correlation coefficients to each prototype
res <- compare.proto.cor(gene.cor=probecor$cor, proto.cor=protocor$cor,
  nn=probecor$cor.n, p.adjust.m="fdr")
head(res)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ correlation }
