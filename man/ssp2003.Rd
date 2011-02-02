\name{ssp2003}
\alias{ssp2003}
\alias{ssp2003.scale}
\alias{ssp2003.robust}
\docType{data}
\title{
SSP2003 classifier for identification of breast cancer molecular subtypes (Sorlie et al 2003)
}
\description{
List of parameters defining the SSP2003 classifier for identification of breast cancer molecular subtypes (Sorlie et al 2003). 
}
\usage{
data(ssp2003)
data(ssp2003.scale)
data(ssp2003.robust)
}
\format{
List of parameters for SSP2003:
\describe{
\item{\code{centroids} }{Gene expression centroids for each subtype.}
\item{\code{centroids.map} }{Mapping for centroids.}
\item{\code{method.cor} }{Method of correlation used to compute distance to the centroids.}
\item{\code{method.centroids} }{Method used to compute the centroids.}
\item{\code{std} }{Method of standardization for gene expressions.}
\item{\code{mins} }{Minimum number of samples within each cluster allowed during the fitting of the model.}
}
}
\details{ 
Three versions of the model are provided, each of ones differs by the gene expressions standardization method since it has an important impact on the subtype classification:
\describe{
\item{\code{ssp2003} }{Use of the official centroids without scaling of the gene expressions.}
\item{\code{ssp2003.scale} }{Use of the official centroids with traditional scaling of the gene expressions (see \code{\link[base]{scale}}).}
\item{\code{ssp2003.robust} }{Use of the official centroids with robust scaling of the gene expressions (see \code{\link[genefu]{rescale}}).}
}
The model \code{ssp2003.robust} has been shown to reach the best concordance with the traditional clinical parameters (ER IHC, HER2 IHC/FISH and histological grade). However the use of this model is recommended only when the dataset is representative of a global population of breast cancer patients (no sampling bias, the 5 subtypes should be present).
}
\source{
\url{http://www.pnas.org/content/100/14/8418}
}
\references{
T. Sorlie and R. Tibshirani and J. Parker and T. Hastie and J. S. Marron and A. Nobel and S. Deng and H. Johnsen and R. Pesich and S. Geister and J. Demeter and C. Perou and P. E. Lonning and P. O. Brown and A. L. Borresen-Dale and D. Botstein (2003) "Repeated Observation of Breast Tumor Subtypes in Independent Gene Expression Data Sets", \emph{Proceedings of the National Academy of Sciences}, \bold{1}(14):8418--8423
}
\examples{
data(ssp2003)
str(ssp2003)
data(ssp2003.robust)
str(ssp2003.robust)
}
\keyword{datasets}
