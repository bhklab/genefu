\name{spearmanCI}
\alias{spearmanCI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the confidence interval for the Spearman correelation coefficient
}
\description{
This function enables to compute the confidence interval for the Spearman correelation coefficient using the Fischer Z transformation
}
\usage{
spearmanCI(x, n, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Spearman correlation coefficient rho
}
  \item{n}{
the sample size used to compute the Spearman rho
}
  \item{alpha}{
alpha level for confidence interval
}
}
%%\details{
%%}
\value{
a vector containing the lower, upper values for the confidence interval and p-value for Spearman rho
}
%%\references{
%%}
\author{
Benjamin Haibe-Kains
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%%}
\examples{
spearmanCI(x=0.2, n=100, alpha=0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
