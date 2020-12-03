\name{collapseIDs}
\alias{collapseIDs}
\title{Utility function to collapse IDs}
\description{
 Utility function called within the claudinLow classifier
}
\usage{
 collapseIDs(x,allids=row.names(x),method="mean")
}
\arguments{
  \item{x}{Matrix of numbers}
  \item{allids}{Defaults to rownames of matrix}
  \item{method}{Default method is "mean"}
}
\references{
  \code{citation("claudinLow")}
}
\seealso{ \code{\link[genefu]{claudinLow}} }


