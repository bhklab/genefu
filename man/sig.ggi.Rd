\name{sig.ggi}
\alias{sig.ggi}
\docType{data}
\title{
Gene expression Grade Index (GGI) as published in Sotiriou et al. 2006
}
\description{
List of 128 affymetrix hgu133a probesets representing 97 unique genes included in the GGI signature. The "EntrezGene.ID" column allows for mapping and "grade" defines the up-regulation of the expressions either in histological grade 1 or 3.
}
\usage{data(sig.ggi)}
\format{
  \code{\link[genefu]{sig.ggi}} is a matrix with 9 columns containing the annotations and information related to the signature itself.
}
%%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%%}
\source{
\url{http://jnci.oxfordjournals.org/cgi/content/full/98/4/262/DC1}
}
\references{
Sotiriou C, Wirapati P, Loi S, Harris A, Bergh J, Smeds J, Farmer P, Praz V, Haibe-Kains B, Lallemand F, Buyse M, Piccart MJ and Delorenzi M (2006) "Gene expression profiling in breast cancer: Understanding the molecular basis of histologic grade to improve prognosis", \emph{Journal of National Cancer Institute}, \bold{98}:262--272
}
\examples{
data(sig.ggi)
head(sig.ggi)
}
\keyword{data}
