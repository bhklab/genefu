\name{nkis}
\alias{nkis}
\alias{data.nkis}
\alias{annot.nkis}
\alias{demo.nkis}
\docType{data}
\title{
Gene expression, annotations and clinical data from van de Vijver et al. 2002
}
\description{
This dataset contains (part of) the gene expression, annotations and clinical data as published in van de Vijver et al. 2002.
}
\usage{data(nkis)}
\format{
\code{\link[genefu]{nkis}} is a dataset containing three matrices:
 \describe{
    \item{data.nkis}{Matrix containing gene expressions as measured by Agilent technology (dual-channel, oligonucleotides}
    \item{annot.nkis}{Matrix containing annotations of Agilent microarray platform}
    \item{demo.nkis}{Clinical information of the breast cancer patients whose tumors were hybridized}
}
}
\details{
This dataset represent only partially the one published by van  de Vijver et al. in 2008. Indeed, only part of the patients (150) and gene expressions (922) are contained in \code{\link[genefu]{data.nkis}}.
}
\source{
\url{http://www.nature.com/nature/journal/v415/n6871/full/415530a.html}
}
\references{
M. J. van de Vijver and Y. D. He and L. van't Veer and H. Dai and A. M. Hart and D. W. Voskuil and G. J. Schreiber and J. L. Peterse and C. Roberts and M. J. Marton and M. Parrish and D. Atsma and A. Witteveen and A. Glas and L. Delahaye and T. van der Velde and H. Bartelink and S. Rodenhuis and E. T. Rutgers and S. H. Friend and R. Bernards (2002) "A Gene Expression Signature as a Predictor of Survival in Breast Cancer", \emph{New England Journal of Medicine}, \bold{347}(25):1999--2009
}
\examples{
data(nkis)
}
\keyword{data}
