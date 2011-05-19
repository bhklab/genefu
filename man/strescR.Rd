\name{strescR}
\alias{strescR}
\title{Utility function to escape LaTeX special characters present in a string}
\description{
 This function returns a vector of strings in which LaTeX special characters are
 escaped, this was useful in conjunction with xtable.
}
\usage{
strescR(strings)
}
\arguments{
  \item{strings}{A vector of strings to deal with.}
}
\value{
  Returns a vector of strings with escaped characters within each string.
}
\references{
  \code{citation("seqinr")}
}
\author{J.R. Lobry}
\seealso{ \code{\link[seqinr]{stresc}} }
\examples{
  strescR("MISC_RNA")
  strescR(c("BB_0001","BB_0002"))
}


