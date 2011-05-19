\name{read.m.file}
\alias{read.m.file}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to read a 'csv' file containing gene lists (aka gene signatures)
}
\description{
This function allows for reading a 'csv' file containing gene signatures. Each gene signature is composed of at least four columns: "gene.list" is the name of the signature on the first line and empty fields below, "probes" are the probe names, "EntrezGene.ID" are the EntrezGene IDs and "coefficient" are the coefficients of each probe.
}
\usage{
read.m.file(file, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{ Filename of the 'csv' file. }
  \item{\dots}{ Additional parameters for \code{\link[utils]{read.csv}} function. }
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
%%  ~Describe the value returned
List of gene signatures.
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
\code{\link[genefu]{mod1}}, \code{\link[genefu]{mod2}}, 'extdata/desmedt2008_genemodules.csv', 'extdata/haibekains2009_sig_genius.csv'
}
\examples{
## read the seven gene modules as published in Desmedt et al 2008
genemods <- read.m.file(system.file("extdata/desmedt2008_genemodules.csv",
  package = "genefu"))
str(genemods, max.level=1)
## read the three subtype signtaures from GENIUS
geniusm <- read.m.file(system.file("extdata/haibekains2009_sig_genius.csv",
  package = "genefu"))
str(geniusm, max.level=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }

