\name{weighted.ttest}
\alias{weighted.ttest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the weighted mean and weighted variance of 'x'
}
\description{
This function allows for computing the weighted mean and weighted variance of a vector of continuous values.
}
\usage{
weighted.ttest(x, w1, w2, alternative=c("two.sided", "less", "greater"), na.rm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object containing the values whose weighted mean is to be computed.
}
  \item{w1}{
a numerical vector of weights of the same length as \code{x} giving the weights to use for elements of \code{x} in the first class.
}
  \item{w2}{
a numerical vector of weights of the same length as \code{x} giving the weights to use for elements of \code{x} in the second class.
}
  \item{alternative}{
a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".  You can specify just the initial letter.
}
  \item{na.rm}{
\code{TRUE} if missing values should be removed, \code{FALSE} otherwise.
}
}
\details{
The weights \code{w1} and \code{w2} should represent the likelihood for each observation stored in \code{x} to belong to the first and second class, respectively. Therefore the values contained in \code{w1} and \code{w2} should lay in [0,1] and 0 <= (w1[i] + w2[i]) <= 1 for i in \{0,1,...,n\} where n is the length of x.

The Welch's version of the t test is implemented in this function, therefore assuming unequal sample size and unequal variance. The sample size of the first and second class are calculated as the \code{sum(w1)} and \code{sum(w2)}, respectively.
}
\value{
A numeric vector of six values that are the difference between the two weighted means, the value of the t statistic, the sample size of class 1, the sample size of class 2, the degree of freedom and the corresponding p-value.
}
\references{
\url{http://en.wikipedia.org/wiki/T_test}
\url{http://www.nicebread.de/blog/files/fc02e1635792cb0f2b3cbd1f7e6c580b-10.php}
}
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
## random generation of 50 normally distributed values for each of the two classes
xx <- c(rnorm(50), rnorm(50)+1)
## fuzzy membership to class 1
ww1 <- runif(50) + 0.3
ww1[ww1 > 1] <- 1
ww1 <- c(ww1, 1 - ww1)
## fuzzy membership to class 2
ww2 <- 1 - ww1
## Welch's t test weighted by fuzzy membership to class 1 and 2 
weighted.ttest(x=xx, w1=ww1, w2=ww2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
