\name{power.cor}
\alias{power.cor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for sample size calculation for correlation coefficients
}
\description{
This function enables to compute the sample size requirements for estimating pearson, kendall and spearman correlations
}
\usage{
power.cor(rho, w, alpha = 0.05, method = c("pearson", "kendall", "spearman"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rho}{
Correaltion coefficients rho (Pearson, Kendall or Spearman)
}
  \item{w}{
a numerical vector of weights of the same length as \code{x} giving the weights to use for elements of \code{x} in the first class.
}
  \item{alpha}{
alpha level
}
  \item{method}{
a character string specifying the method to compute the correlation coefficient, must be one of "pearson" (default), "kendall" or "spearman".  You can specify just the initial letter.
}
}
%%\details{
%%}
\value{
sample size requirement
}
\references{
Bonett, D. G., and Wright, T. A. (2000). Sample size requirements for estimating pearson, kendall and spearman correlations. Psychometrika, 65(1), 23-28. doi:10.1007/BF02294183
}
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
power.cor(rho=0.5, w=0.1, alpha=0.05, method="spearman")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
