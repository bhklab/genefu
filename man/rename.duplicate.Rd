\name{rename.duplicate}
\alias{rename.duplicate}
\title{
Function to rename duplicated strings.
}
\description{
This function renames duplicated strings by adding their number of occurrences at the end.
}
\usage{
rename.duplicate(x, sep = "_", verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ vector of strings. }
  \item{sep}{ a character to be the separator between the number added at the end and the string itself. }
  \item{verbose}{ \code{TRUE} to print informative messages, \code{FALSE} otherwise. }
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
  \item{new.x }{new strings (without duplicates). }
  \item{duplicated.x }{ strings which were originally duplicated. }
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
nn <- sample(letters[1:10], 30, replace=TRUE)
table(nn)
rename.duplicate(x=nn, verbose=TRUE)
}
\keyword{ character }
