\name{map.datasets}
\alias{map.datasets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to map a list of datasets through EntrezGene IDs in order to get the union of the genes
}
\description{
This function maps a list of datasets through EntrezGene IDs in order to get the union of the genes. 
}
\usage{
map.datasets(datas, annots, do.mapping = FALSE, mapping, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datas}{
List of matrices of gene expressions with samples in rows and probes in columns, dimnames being properly defined.
}
  \item{annots}{
List of matrices of annotations with at least one column named "EntrezGene.ID", dimnames being properly defined.
}
  \item{do.mapping}{
\code{TRUE} if the mapping through Entrez Gene ids must be performed (in case of ambiguities, the most variant probe is kept for each gene), \code{FALSE} otherwise.
}
  \item{mapping}{
Matrix with columns "EntrezGene.ID" and "probe.x" used to force the mapping such that the probes of platform x are not selected based on their variance.
}
  \item{verbose}{
\code{TRUE} to print informative messages, \code{FALSE} otherwise.
}
}
\details{
In case of several probes representing the same EntrezGene ID, the most variant is selected if \code{mapping} is not specified. When a EntrezGene ID does not exist in a specific dataset, NA values are introduced.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{datas }{List of datasets (gene expression matrices)}
  \item{annots }{List of annotations (annotation matrices)}
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

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
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
str(datas.mapped, max.level=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ mapping }

