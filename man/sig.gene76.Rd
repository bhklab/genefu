\name{sig.gene76}
\alias{sig.gene76}
\docType{data}
\title{
Signature used to compute the Relapse Score (GENE76) as published in Wang et al. 2005
}
\description{
List of 76 affymetrix hgu133a probesets representing 60 unique genes included in the GENE76 signature. The EntrezGene.ID allows for mapping and the coefficient allows for signature computation.
}
\usage{data(sig.gene76)}
\format{
  \code{\link[genefu]{sig.gene76}} is a matrix with 10 columns containing the annotations and information related to the signature itself.
}
%%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%%}
\source{
\url{http://www.thelancet.com/journals/lancet/article/PIIS0140-6736(05)17947-1/abstract}
}
\references{
Y. Wang and J. G. Klijn and Y. Zhang and A. M. Sieuwerts and M. P. Look and F. Yang and D. Talantov and M. Timmermans and M. E. Meijer-van Gelder and J. Yu and T. Jatkoe and E. M. Berns and D. Atkins and J. A. Foekens (2005) "Gene-Expression Profiles to Predict Distant Metastasis of Lymph-Node-Negative Primary Breast Cancer", \emph{Lancet}, \bold{365}(9460):671--679.
}
\examples{
data(sig.gene76)
head(sig.gene76)
}
\keyword{data}
