\name{sig.endoPredict}
\alias{sig.endoPredict}
\docType{data}
\title{
Signature used to compute the endoPredict signature as published by Filipits et al 2011
}
\description{
List of 11 genes included in the endoPredict signature. The EntrezGene.ID allows for mapping and the mapping to affy probes is already provided.
}
\usage{data(sig.endoPredict)}
\format{
  \code{\link[genefu]{sig.endoPredict}} is a matrix with 5 columns containing the annotations and information related to the signature itself (including a mapping to Affymetrix HGU platform).
}
%%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%%}
%%\source{
%%}
\references{
Filipits, M., Rudas, M., Jakesz, R., Dubsky, P., Fitzal, F., Singer, C. F., et al. (2011). "A new molecular predictor of distant recurrence in ER-positive, HER2-negative breast cancer adds independent information to conventional clinical risk factors." \emph{Clinical Cancer Research}, \bold{17}(18):6012-–6020.
}
\examples{
data(sig.endoPredict)
head(sig.endoPredict)
}
\keyword{data}
