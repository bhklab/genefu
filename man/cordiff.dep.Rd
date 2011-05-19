\name{cordiff.dep}
\alias{cordiff.dep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to estimate whether two dependent correlations differ
}
\description{
This function tests for statistical differences between two dependent correlations using the formula provided on page 56 of Cohen & Cohen (1983). The function returns a t-value, the DF and the p-value.
}
\usage{
cordiff.dep(r.x1y, r.x2y, r.x1x2, n,
  alternative = c("two.sided", "less", "greater"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r.x1y}{
The correlation between x1 and y where y is typically your outcome variable.
}
  \item{r.x2y}{
The correlation between x2 and y where y is typically your outcome variable.
}
  \item{r.x1x2}{
The correlation between x1 and x2 (the correlation between your two predictors).
}
  \item{n}{
The sample size.
}
  \item{alternative}{
A character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
}
}
\details{
This function is inspired from the \code{\link[multilevel]{cordif.dep}}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
Vector of three values: t statistics, degree of freedom, and p-value.
}
\references{
Cohen, J. & Cohen, P. (1983) "Applied multiple regression/correlation analysis for the behavioral sciences (2nd Ed.)" \emph{Hillsdale, nJ: Lawrence Erlbaum Associates}.
}
\author{
Benjamin Haibe-Kains
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[stats]{cor}}, \code{\link[stats]{t.test}}, \code{\link[genefu]{compare.proto.cor}}
}
\examples{
## load VDX dataset
data(vdx)
## retrieve ESR1, AURKA and MKI67 gene expressions
x1 <- data.vdx[ ,"208079_s_at"]
x2 <- data.vdx[ ,"205225_at"]
y <- data.vdx[ ,"212022_s_at"]
## is MKI67 significantly more correlated to AURKA than ESR1?
cc.ix <- complete.cases(x1, x2, y)
cordiff.dep(r.x1y=abs(cor(x=x1[cc.ix], y=y[cc.ix], use="everything",
  method="pearson")), r.x2y=abs(cor(x=x2[cc.ix], y=y[cc.ix],
  use="everything", method="pearson")), r.x1x2=abs(cor(x=x1[cc.ix],
  y=x2[cc.ix], use="everything", method="pearson")), n=sum(cc.ix),
  alternative="greater")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
