\name{boxplotplus2}
\alias{boxplotplus2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Box plot of group of values with corresponding jittered points
}
\description{
This function allows for display a boxplot with jittered points.
}
\usage{
boxplotplus2(x, .jit = 0.25, .las = 1, .ylim, box.col = "lightgrey",
  pt.col = "blue", pt.cex = 0.5, pt.pch = 16, med.line = FALSE,
  med.col = "goldenrod", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
\code{x} could be a list of group values or a matrix (each group is a row).
}
  \item{.jit}{
Amount of jittering noise.
}
  \item{.las}{
Numeric in {0,1,2,3}; the style of axis labels.
}
  \item{.ylim}{
Range for y axis.
}
  \item{box.col}{
Color for boxes.
}
  \item{pt.col}{
Color for groups (jittered points).
}
  \item{pt.cex}{
A numerical value giving the amount by which plotting jittered points should be magnified relative to the default.
}
  \item{pt.pch}{
Either an integer specifying a symbol or a single character to be used as the default in plotting jittered points.  See \code{\link[graphics]{points}} for possible values and their interpretation.
}
  \item{med.line}{
\code{TRUE} if a line should link the median of each group, \code{FALSE} otherwise.
}
  \item{med.col}{
Color of \code{med.line}.
}
  \item{\dots}{
Additional parameters for \code{\link[graphics]{boxplot}} function.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
%%  ~Describe the value returned
Number of samples in each group.
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Christos Hatzis
}
\note{
2.21.2006 - Christos Hatzis, Nuvera Biosciences
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[graphics]{boxplot}}, \code{\link[base]{jitter}}
}
\examples{
dd <- list("G1"=runif(20), "G2"=rexp(30) * -1.1, "G3"=rnorm(15) * 1.3)
boxplotplus2(x=dd, .las=3, .jit=0.75, .ylim=c(-3,3), pt.cex=0.75,
  pt.col=c(rep("darkred", 20), rep("darkgreen", 30), rep("darkblue", 15)),
  pt.pch=c(0, 9, 17))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
