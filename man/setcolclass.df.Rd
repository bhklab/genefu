\name{setcolclass.df}
\alias{setcolclass.df}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to set the class of columns in a data.frame
}
\description{
This function enables to set the class of each culumn in a data.frame
}
\usage{
setcolclass.df(df, colclass, factor.levels)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
data.frame for which columns' class need to be updated
}
  \item{colclass}{
class for each column of the data.frame
}
  \item{factor.levels}{
list of levels for each factor
}
}
%%\details{
%%}
\value{
A data.frame with columns' class and levels properly set
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
tt <- data.frame(matrix(NA, nrow=3, ncol=3, dimnames=list(1:3, paste("column", 1:3, sep="."))), 
stringsAsFactors=FALSE)
tt <- setcolclass.df(df=tt, colclass=c("numeric", "factor", "character"), 
factor.levels=list(NULL, c("F1", "F2", "F3"), NULL))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{  }
