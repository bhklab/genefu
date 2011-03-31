\name{pam50}
\alias{pam50}
\alias{pam50.scale}
\alias{pam50.robust}
\docType{data}
\title{
PAM50 classifier for identification of breast cancer molecular subtypes (Parker et al 2009)
}
\description{
List of parameters defining the PAM50 classifier for identification of breast cancer molecular subtypes (Parker et al 2009). 
}
\usage{
data(pam50)
data(pam50.scale)
data(pam50.robust)
}
\format{
List of parameters for PAM50:
\describe{
\item{\code{centroids} }{Gene expression centroids for each subtype.}
\item{\code{centroids.map} }{Mapping for centroids.}
\item{\code{method.cor} }{Method of correlation used to compute distance to the centroids.}
\item{\code{method.centroids} }{Method used to compute the centroids.}
\item{\code{std} }{Method of standardization for gene expressions ("none", "scale" or "robust").}
\item{\code{mins} }{Minimum number of samples within each cluster allowed during the fitting of the model.}
}
}
\details{ 
Three versions of the model are provided, each of ones differs by the gene expressions standardization method since it has an important impact on the subtype classification:
\describe{
\item{\code{pam50} }{Use of the official centroids without scaling of the gene expressions.}
\item{\code{pam50.scale} }{Use of the official centroids with traditional scaling of the gene expressions (see \code{\link[base]{scale}}).}
\item{\code{pam50.robust} }{Use of the official centroids with robust scaling of the gene expressions (see \code{\link[genefu]{rescale}}).}
}
The model \code{pam50.robust} has been shown to reach the best concordance with the traditional clinical parameters (ER IHC, HER2 IHC/FISH and histological grade). However the use of this model is recommended only when the dataset is representative of a global population of breast cancer patients (no sampling bias, the 5 subtypes should be present).
}
\source{
\url{http://jco.ascopubs.org/cgi/content/short/JCO.2008.18.1370v1}
}
\references{
Parker, Joel S. and Mullins, Michael and Cheang, Maggie C.U. and Leung, Samuel and Voduc, David and Vickery, Tammi and Davies, Sherri and Fauron, Christiane and He, Xiaping and Hu, Zhiyuan and Quackenbush, John F. and Stijleman, Inge J. and Palazzo, Juan and Marron, J.S. and Nobel, Andrew B. and Mardis, Elaine and Nielsen, Torsten O. and Ellis, Matthew J. and Perou, Charles M. and Bernard, Philip S. (2009) "Supervised Risk Predictor of Breast Cancer Based on Intrinsic Subtypes", \emph{Journal of Clinical Oncology}, \bold{27}(8):1160--1167
}
\examples{
data(pam50)
str(pam50)
data(pam50.robust)
str(pam50.robust)
}
\keyword{datasets}
