\name{weighted.meanvar}
\alias{weighted.meanvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the weighted mean and weighted variance of 'x'
}
\description{
This function allows for computing the weighted mean and weighted variance of a vector of continuous values.
}
\usage{
weighted.meanvar(x, w, na.rm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object containing the values whose weighted mean is to be computed.
}
  \item{w}{
a numerical vector of weights of the same length as \code{x} giving the weights to use for elements of \code{x}.
}
  \item{na.rm}{
\code{TRUE} if missing values should be removed, \code{FALSE} otherwise.
}
}
\details{
If \code{w} is missing then all elements of \code{x} are given the same weight, otherwise the weights coerced to numeric by \code{as.numeric}. On the contrary of \link[stats]{weighted.mean} the weights are NOT normalized to sum to one. If the sum of the weights is zero or infinite, NAs will be returned.
}
\value{
A numeric vector of two values that are the weighted mean and weighted variance respectively.
}
\references{
\url{http://en.wikipedia.org/wiki/Weighted_variance#Weighted_sample_variance}}
\author{
Benjamin Haibe-Kains
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link[stats]{weighted.mean}
}
\examples{
set.seed(54321)
weighted.meanvar(x=rnorm(100) + 10, w=runif(100))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univar }
