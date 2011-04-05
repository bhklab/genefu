\name{stab.fs.ranking}
\alias{stab.fs.ranking}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to quantify stability of feature ranking.
}
\description{
This function computes several indexes to quantify feature ranking stability for several number of selected features. This is usually estimated through perturbation of the original dataset by generating multiple sets of selected features.
}
\usage{
stab.fs.ranking(fsets, sizes, N, method = c("kuncheva", "davis"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fsets}{ list or matrix of sets of selected features (in rows), each ranking must have the same size }
  \item{sizes}{ Number of top-ranked features for which the stability index must be computed
}
  \item{N}{ total number of features on which feature selection is performed }
  \item{method}{ stability index (see details section) }
  \item{\dots}{ additional parameters passed to stability index (\code{penalty} that is a numeric for Davis' stability index, see details section) }
}
\details{
Stability indices may use different parameters. In this version only the Davis index requires an additional parameter that is \code{penalty}, a numeric value used as penalty term.

Kuncheva index (\code{kuncheva}) lays in [-1, 1], An index of -1 means no intersection between sets of selected features, +1 means that all the same features are always selected and 0 is the expected stability of a random feature selection.

Davis index (\code{davis}) lays in [0,1], With a pnalty term equal to 0, an index of 0 means no intersection between sets of selected features and +1 means that all the same features are always selected. A penalty of 1 is usually used so that a feature selection performed with no or all features has a Davis stability index equals to 0. None estimate of the expected Davis stability index of a random feature selection was published.
}
\value{ A vector of numeric that are stability indices for each size of the sets of selected features given the rankings
}
\references{
Davis CA, Gerick F, Hintermair V, Friedel CC, Fundel K, Kuffner R, Zimmer R (2006) "Reliable gene signatures for microarray classification: assessment of stability and performance", \emph{Bioinformatics}, \bold{22}(19):356-2363.

Kuncheva LI (2007) "A stability index for feature selection", \emph{AIAP'07: Proceedings of the 25th conference on Proceedings of the 25th IASTED International Multi-Conference}, pages 390--395.
}
\author{
Benjamin Haibe-Kains
}
%%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[genefu]{stab.fs}}
}
\examples{
## 100 random selection of 50 features from a set of 10,000 features
fsets <- lapply(as.list(1:100), function(x, size=50, N=10000) {
  return(sample(1:N, size, replace=FALSE))} )
names(fsets) <- paste("fsel", 1:length(fsets), sep=".")

## Kuncheva index
stab.fs.ranking(fsets=fsets, sizes=c(1, 10, 20, 30, 40, 50),
  N=10000, method="kuncheva")
## close to 0 as expected for a random feature selection

## Davis index
stab.fs.ranking(fsets=fsets, sizes=c(1, 10, 20, 30, 40, 50),
  N=10000, method="davis", penalty=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ feature selection }
\keyword{ stability }% __ONLY ONE__ keyword per line
