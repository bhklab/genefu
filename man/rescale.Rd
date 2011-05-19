\name{rescale}
\alias{rescale}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to rescale values based on quantiles
}
\description{
This function rescales values \code{x} based on quantiles specified by the user such that x' = (x - q1) / (q2 - q1) where q is the specified quantile, q1 = q / 2, q2 = 1 - q/2) and x' are the new rescaled values.
}
\usage{
rescale(x, na.rm = FALSE, q = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{na.rm}{
\code{TRUE} if missing values should be removed, \code{FALSE} otherwise.
}
  \item{q}{
Quantile (must lie in [0,1]).
}
}
\details{
In order to rescale gene expressions, \code{q = 0.05} yielded comparable scales in numerous breast cancer microarray datasets (data not shown).The rational behind this is that, in general, 'extreme cases' (e.g. low and high proliferation, high and low expression of ESR1, ...) are often present in microarray datasets, making the estimation of 'extreme' quantiles quite stable. This is specially true for genes exhibiting some multi-modality like ESR1 or ERBB2.
}
\value{
%%  ~Describe the value returned
Vector of rescaled values with two attributes \code{q1} and \code{q1} containing the values of the lower and the upper quantiles respectively.
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
\code{\link[base]{scale}}
}
\examples{
## load VDX dataset
data(vdx)
## load NKI dataset
data(nki)
## example of rescaling for ESR1 expression
par(mfrow=c(2,2))
hist(data.vdx[ ,"205225_at"], xlab="205225_at", breaks=20,
  main="ESR1 in VDX")
hist(data.nki[ ,"NM_000125"], xlab="NM_000125", breaks=20,
  main="ESR1 in NKI")
hist((rescale(x=data.vdx[ ,"205225_at"], q=0.05) - 0.5) * 2,
  xlab="205225_at", breaks=20, main="ESR1 in VDX\nrescaled")
hist((rescale(x=data.nki[ ,"NM_000125"], q=0.05) - 0.5) * 2,
  xlab="NM_000125", breaks=20, main="ESR1 in NKI\nrescaled")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
