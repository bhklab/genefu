\name{sig.oncotypedx}
\alias{sig.oncotypedx}
\docType{data}
\title{
Signature used to compute the OncotypeDX signature as published by Paik et al 2004
}
\description{
List of 21 genes included in the OncotypeDX signature. The EntrezGene.ID allows for mapping and the mapping to affy probes is already provided.
}
\usage{data(sig.oncotypedx)}
\format{
  \code{\link[genefu]{sig.oncotypedx}} is a matrix with 5 columns containing the annotations and information related to the signature itself (including a mapping to Affymetrix HGU platform).
}
%%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%%}
%%\source{
%%}
\references{
S. Paik, S. Shak, G. Tang, C. Kim, J. Bakker, M. Cronin, F. L. Baehner, M. G. Walker, D. Watson, T. Park, W. Hiller, E. R. Fisher, D. L. Wickerham, J. Bryant, and N. Wolmark (2004) "A Multigene Assay to Predict Recurrence of Tamoxifen-Treated, Node-Negative Breast Cancer", \emph{New England Journal of Medicine}, \bold{351}(27):2817--2826.
}
\examples{
data(sig.oncotypedx)
head(sig.oncotypedx)
}
\keyword{data}
