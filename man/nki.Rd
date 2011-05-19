\name{nki}
\alias{nki}
\alias{data.nki}
\alias{annot.nki}
\alias{demo.nki}
\docType{data}
\title{
Gene expression, annotations and clinical data from van de Vijver et al. 2002
}
\description{
This dataset contains (part of) the gene expression, annotations and clinical data as published in van de Vijver et al. 2002.
}
\usage{data(nki)}
\format{
\code{\link[genefu]{nki}} is a dataset containing three matrices:
 \describe{
    \item{data.nki}{Matrix containing gene expressions as measured by Agilent technology (dual-channel, oligonucleotides}
    \item{annot.nki}{Matrix containing annotations of Agilent microarray platform}
    \item{demo.nki}{Clinical information of the breast cancer patients whose tumors were hybridized}
}
}
\details{
This dataset represent only partially the one published by van  de Vijver et al. in 2008. Indeed, only part of the patients (150) and gene expressions (922) are contained in \code{\link[genefu]{data.nki}}.
}
\source{
\url{http://www.rii.com/publications/2002/vantveer.html}
}
\references{
M. J. van de Vijver and Y. D. He and L. van't Veer and H. Dai and A. M. Hart and D. W. Voskuil and G. J. Schreiber and J. L. Peterse and C. Roberts and M. J. Marton and M. Parrish and D. Atsma and A. Witteveen and A. Glas and L. Delahaye and T. van der Velde and H. Bartelink and S. Rodenhuis and E. T. Rutgers and S. H. Friend and R. Bernards (2002) "A Gene Expression Signature as a Predictor of Survival in Breast Cancer", \emph{New England Journal of Medicine}, \bold{347}(25):1999--2009
}
\examples{
data(nki)
}
\keyword{data}
