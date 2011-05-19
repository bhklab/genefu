\name{tbrm}
\alias{tbrm}
\title{ Function to compute Tukey's Biweight Robust Mean }
\description{
Computation of Tukey's Biweight Robust Mean, a robust average that is unaffected by outliers.
}
\usage{
tbrm(x, C = 9)
}
\arguments{
  \item{x}{ a numeric vector }
  \item{C}{ a constant. \code{C} is preassigned a value of 9 according to the
    Cook reference below but other values are possible. }
}
\details{
  This is a one step computation that follows the Affy whitepaper below see
  page 22. This function is called by \code{\link[dplR]{chron}} to calculate a
  robust mean.  \code{C} determines the point at which outliers are given a
  weight of 0 and therefore do not contribute to the calculation of the mean.
  C=9 sets values roughly +/-6 standard deviations to 0. C=6 is also used in
  tree-ring chronology development. Cook and Kairiukstis (1990) have further
  details.

Retrieved from \code{\link[dplR]{tbrm}}.
}
\value{
  A numeric mean.
}
\references{
Statistical Algorithms Description Document, 2002, Affymetrix. p22.

Cook, E. R. and Kairiukstis, L.A. (1990) \emph{Methods of Dendrochronology: Applications in the Environmental Sciences.} ISBN-13: 978-0792305866.

Mosteller, F. and Tukey, J. W. (1977) \emph{Data Analysis and Regression: a second course in statistics.} Addison-Wesley. ISBN-13: 978-0201048544.
}
\author{ Andy Bunn }
\seealso{ \code{\link[dplR]{chron}} }
\examples{
tbrm(rnorm(100))

}
\keyword{ misc }
