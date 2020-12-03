\name{write.m.file}
\alias{write.m.file}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to write a 'csv' file containing gene lists (aka gene signatures)
}
\description{
This function allows for writing a 'csv' file containing gene signatures. Each gene signature is composed of at least four columns: "gene.list" is the name of the signature on the first line and empty fields below, "probes" are the probe names, "EntrezGene.ID" are the EntrezGene IDs and "coefficient" are the coefficients of each probe.
}
\usage{
write.m.file(obj, file, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
List of gene signatures.
}
  \item{file}{
Filename of the 'csv' file.
}
  \item{\dots}{
Additional parameters for \code{\link[utils]{read.csv}} function.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
None.
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
## load gene modules published by Demsedt et al 2009
data(mod1)
## write these gene modules in a 'csv' file
\dontrun{write.m.file(obj=mod1, file="desmedt2009_genemodules.csv")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
