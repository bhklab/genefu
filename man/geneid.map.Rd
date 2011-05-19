\name{geneid.map}
\alias{geneid.map}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to find the common genes between two datasets or a dataset and a gene list
}
\description{
This function allows for fast mapping between two datasets or a dataset and a gene list. The mapping process is performed using Entrez Gene id as reference. In case of ambiguities (several probes representing the same gene), the most variant probe is selected.
}
\usage{
geneid.map(geneid1, data1, geneid2, data2, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geneid1}{
first vector of Entrez Gene ids. The name of the vector cells must be the name of the probes in the dataset \code{data1}.
}
  \item{data1}{
First dataset with samples in rows and probes in columns. The dimnames must be properly defined.
}
  \item{geneid2}{
Second vector of Entrez Gene ids. The name of the vector cells must be the name of the probes in the dataset \code{data1} if it is not missing, proper names must be assigned otherwise.
}
  \item{data2}{
First dataset with samples in rows and probes in columns. The dimnames must be properly defined. It may be missing.
}
  \item{verbose}{
\code{TRUE} to print informative messages, \code{FALSE} otherwise.
}
}
%%\details{
%%}
\value{
  \item{geneid1 }{Mapped gene list from \code{geneid1}.}
  \item{data1 }{Mapped dataset from \code{data1}.}
  \item{geneid2 }{Mapped gene list from \code{geneid2}.}
  \item{data2 }{Mapped dataset from \code{data2}.}
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Benjamin Haibe-Kains
}
\note{
It is mandatory that the names of \code{geneid1} and \code{geneid2} must be the probe names of the microarray platform.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
## load NKI data
data(nki)
nki.gid <- annot.nki[ ,"EntrezGene.ID"]
names(nki.gid) <- dimnames(annot.nki)[[1]]
## load GGI signature
data(sig.ggi)
ggi.gid <- sig.ggi[ ,"EntrezGene.ID"]
names(ggi.gid) <- as.character(sig.ggi[ ,"probe"])
## mapping through Entrez Gene ids of NKI and GGI signature
res <- geneid.map(geneid1=nki.gid, data1=data.nki,
  geneid2=ggi.gid, verbose=FALSE)
str(res)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ mapping }
