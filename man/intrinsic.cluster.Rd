\name{intrinsic.cluster}
\alias{intrinsic.cluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to fit a Single Sample Predictor (SSP) as in Perou, Sorlie, Hu, and Parker publications
}
\description{
This function fits the Single Sample Predictor (SSP) as published in Sorlie et al 2003, Hu et al 2006 and Parker et al 2009. This model is actually a nearest centroid classifier where the centroids representing the breast cancer molecular subtypes are identified through hierarchical clustering using an "intrinsic gene list".
}
\usage{
intrinsic.cluster(data, annot, do.mapping = FALSE, mapping,
  std = c("none", "scale", "robust"), rescale.q = 0.05, intrinsicg,
  number.cluster = 3, mins = 5, method.cor = c("spearman", "pearson"),
  method.centroids = c("mean", "median", "tukey"), filen, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Matrix of gene expressions with samples in rows and probes in columns, dimnames being properly defined.
}
  \item{annot}{
Matrix of annotations with at least one column named "EntrezGene.ID", dimnames being properly defined.
}
  \item{do.mapping}{
\code{TRUE} if the mapping through Entrez Gene ids must be performed (in case of ambiguities, the most variant probe is kept for each gene), \code{FALSE} otherwise.
}
  \item{mapping}{
Matrix with columns "EntrezGene.ID" and "probe" used to force the mapping such that the probes are not selected based on their variance.
}
  \item{std}{
Standardization of gene expressions: \code{scale} for traditional standardization based on mean and standard deviation, \code{robust} for standardization based on the 0.025 and 0.975 quantiles, \code{none} to keep gene expressions unchanged.
}
  \item{rescale.q}{
	Proportion of expected outliers for (\code{robust}) rescaling the gene expressions.
}
  \item{intrinsicg}{
Intrinsic gene lists. May be specified by the user as a matrix wit hat least 2 columns named probe and EntrezGene.ID for the probe names and the corresponding Entrez Gene ids. The intrinsic gene lists published by Sorlie et al. 2003, Hu et al. 2006 and Parker et al. 2009 are stored in \code{ssp2003}, \code{ssp2006} and \code{pam50} respectively.
}
  \item{number.cluster}{
The number of main clusters to be identified by hierarchical clustering.
}
  \item{mins}{
The minimum number of samples to be in a main cluster.
}
  \item{method.cor}{
Correlation coefficient used to identified the nearest centroid. May be \code{spearman} or \code{pearson}.
}
  \item{method.centroids}{
LMethod to compute a centroid from gene expressions of a cluster of samples: \code{mean}, \code{median} or \code{tukey} (Tukey's Biweight Robust Mean).
}
  \item{filen}{
Name of the csv file where the subtype clustering model must be stored.
}
  \item{verbose}{
\code{TRUE} to print informative messages, \code{FALSE} otherwise.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{model }{Subtype Clustering Model (mixture of three Gaussians), like \code{\link[genefu]{scmgene}}, \code{\link[genefu]{scmod1}} and \code{\link[genefu]{scmod2}} when this function is used on VDX dataset (Wang et al. 2005 and Minn et al. 2007) with the gene modules published in the two references cited below.}
\item{subtype }{Subtypes identified by the SSP. For published intrinsic gene lists, subtypes can be either "Basal", "Her2", "LumA", "LumB" or "Normal".}
\item{subtype.proba }{Probabilities  to belong to each subtype estimated from the correlations to each centroid.}
\item{cor }{Correlation coefficient to each centroid.}
}
\references{
T. Sorlie and R. Tibshirani and J. Parker and T. Hastie and J. S. Marron and A. Nobel and S. Deng and H. Johnsen and R. Pesich and S. Geister and J. Demeter and C. Perou and P. E. Lonning and P. O. Brown and A. L. Borresen-Dale and D. Botstein (2003) "Repeated Observation of Breast Tumor Subtypes in Independent Gene Expression Data Sets", \emph{Proceedings of the National Academy of Sciences}, \bold{1}(14):8418--8423

Hu, Zhiyuan and Fan, Cheng and Oh, Daniel and Marron, JS and He, Xiaping and Qaqish, Bahjat and Livasy, Chad and Carey, Lisa and Reynolds, Evangeline and Dressler, Lynn and Nobel, Andrew and Parker, Joel and Ewend, Matthew and Sawyer, Lynda and Wu, Junyuan and Liu, Yudong and Nanda, Rita and Tretiakova, Maria and Orrico, Alejandra and Dreher, Donna and Palazzo, Juan and Perreard, Laurent and Nelson, Edward and Mone, Mary and Hansen, Heidi and Mullins, Michael and Quackenbush, John and Ellis, Matthew and Olopade, Olufunmilayo and Bernard, Philip and Perou, Charles (2006) "The molecular portraits of breast tumors are conserved across microarray platforms", \emph{BMC Genomics}, \bold{7}(96)

Parker, Joel S. and Mullins, Michael and Cheang, Maggie C.U. and Leung, Samuel and Voduc, David and Vickery, Tammi and Davies, Sherri and Fauron, Christiane and He, Xiaping and Hu, Zhiyuan and Quackenbush, John F. and Stijleman, Inge J. and Palazzo, Juan and Marron, J.S. and Nobel, Andrew B. and Mardis, Elaine and Nielsen, Torsten O. and Ellis, Matthew J. and Perou, Charles M. and Bernard, Philip S. (2009) "Supervised Risk Predictor of Breast Cancer Based on Intrinsic Subtypes", \emph{Journal of Clinical Oncology}, \bold{27}(8):1160--1167
}
\author{
Benjamin Haibe-Kains
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[genefu]{subtype.cluster}}, \code{\link[genefu]{intrinsic.cluster.predict}}, \code{\link[genefu]{ssp2003}}, \code{\link[genefu]{ssp2006}}, \code{\link[genefu]{pam50}}
}
\examples{
## load SSP signature published in Sorlie et al. 2003
data(ssp2003)
## load NKI data
data(nki)
## load VDX data
data(vdx)
ssp2003.nki <- intrinsic.cluster(data=data.nki, annot=annot.nki,
  do.mapping=TRUE, std="robust",
  intrinsicg=ssp2003$centroids.map[ ,c("probe", "EntrezGene.ID")],
  number.cluster=5, mins=5, method.cor="spearman",
  method.centroids="mean", verbose=TRUE)
str(ssp2003.nki, max.level=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ clustering }
