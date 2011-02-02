\name{sig.gene70}
\alias{sig.gene70}
\docType{data}
\title{
Signature used to compute the 70 genes prognosis profile (GENE70) as published by van't Veer et al. 2002
}
\description{
List of 70 agilent probe ids representing 56 unique genes included in the GENE70 signature. The EntrezGene.ID allows for mapping and the "average.good.prognosis.profile" values allows for signature computation.
}
\usage{data(sig.gene70)}
\format{
  \code{\link[genefu]{sig.gene70}} is a matrix with 9 columns containing the annotations and information related to the signature itself.
}
%%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%%}
\source{
\url{http://www.rii.com/publications/2002/vantveer.html}
}
\references{
L. J. van't Veer and H. Dai and M. J. van de Vijver and Y. D. He and A. A. Hart and M. Mao and H. L. Peterse and K. van der Kooy and M. J. Marton and A. T. Witteveen and G. J. Schreiber and R. M. Kerkhiven and C. Roberts and P. S. Linsley and R. Bernards and S. H. Friend (2002) "Gene Expression Profiling Predicts Clinical Outcome of Breast Cancer", \emph{Nature}, \bold{415}:530--536.
}
\examples{
data(sig.gene70)
head(sig.gene70)
}
\keyword{data}
