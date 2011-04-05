\name{intrinsic.cluster.predict}
\alias{intrinsic.cluster.predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to identify breast cancer molecular subtypes using the Single Sample Predictor (SSP)
}
\description{
This function identifies the breast cancer molecular subtypes using a Single Sample Predictor (SSP) fitted by \code{\link[genefu]{intrinsic.cluster}}.
}
\usage{
intrinsic.cluster.predict(sbt.model, data, annot, do.mapping = FALSE,
  mapping, do.prediction.strength = FALSE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sbt.model}{
Subtype Clustering Model as returned by \code{\link[genefu]{intrinsic.cluster}}.
}
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
Matrix with columns "EntrezGene.ID" and "probe" used to force the mapping such that the probes are not selected based on  their variance.
}
  \item{do.prediction.strength}{
\code{TRUE} if the prediction strength must be computed (Tibshirani and Walther 2005), \code{FALSE} otherwise.
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
\item{subtype }{Subtypes identified by the SSP. For published intrinsic gene lists, subtypes can be either "Basal", "Her2", "LumA", "LumB" or "Normal".}
\item{subtype.proba }{Probabilities  to belong to each subtype estimated from the correlations to each centroid.}
\item{cor }{Correlation coefficient to each centroid.}
\item{prediction.strength }{Prediction strength for \code{subtypes}.}
\item{subtype.train }{Classification (similar to \code{subtypes}) computed during fitting of the model for prediction strength.}
\item{centroids.map }{Mapped probes from the intrinsic gene list used to compute the centroids.}
\item{profiles }{Intrinsic gene expression profiles for each sample.}
}
\references{
T. Sorlie and R. Tibshirani and J. Parker and T. Hastie and J. S. Marron and A. Nobel and S. Deng and H. Johnsen and R. Pesich and S. Geister and J. Demeter and C. Perou and P. E. Lonning and P. O. Brown and A. L. Borresen-Dale and D. Botstein (2003) "Repeated Observation of Breast Tumor Subtypes in Independent Gene Expression Data Sets", \emph{Proceedings of the National Academy of Sciences}, \bold{1}(14):8418--8423

Hu, Zhiyuan and Fan, Cheng and Oh, Daniel and Marron, JS and He, Xiaping and Qaqish, Bahjat and Livasy, Chad and Carey, Lisa and Reynolds, Evangeline and Dressler, Lynn and Nobel, Andrew and Parker, Joel and Ewend, Matthew and Sawyer, Lynda and Wu, Junyuan and Liu, Yudong and Nanda, Rita and Tretiakova, Maria and Orrico, Alejandra and Dreher, Donna and Palazzo, Juan and Perreard, Laurent and Nelson, Edward and Mone, Mary and Hansen, Heidi and Mullins, Michael and Quackenbush, John and Ellis, Matthew and Olopade, Olufunmilayo and Bernard, Philip and Perou, Charles (2006) "The molecular portraits of breast tumors are conserved across microarray platforms", \emph{BMC Genomics}, \bold{7}(96)

Parker, Joel S. and Mullins, Michael and Cheang, Maggie C.U. and Leung, Samuel and Voduc, David and Vickery, Tammi and Davies, Sherri and Fauron, Christiane and He, Xiaping and Hu, Zhiyuan and Quackenbush, John F. and Stijleman, Inge J. and Palazzo, Juan and Marron, J.S. and Nobel, Andrew B. and Mardis, Elaine and Nielsen, Torsten O. and Ellis, Matthew J. and Perou, Charles M. and Bernard, Philip S. (2009) "Supervised Risk Predictor of Breast Cancer Based on Intrinsic Subtypes", \emph{Journal of Clinical Oncology}, \bold{27}(8):1160--1167

Tibshirani R and Walther G (2005) "Cluster Validation by Prediction Strength", \emph{Journal of Computational and Graphical Statistics}, \bold{14}(3):511--528
}
\author{
Benjamin Haibe-Kains
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[genefu]{intrinsic.cluster}}, \code{\link[genefu]{ssp2003}}, \code{\link[genefu]{ssp2006}}, \code{\link[genefu]{pam50}}
}
\examples{
## load SSP fitted in Sorlie et al. 2003
data(ssp2003)
## load NKI data
data(nki)
## SSP2003 applied on NKI
ssp2003.nki <- intrinsic.cluster.predict(sbt.model=ssp2003,
  data=data.nki, annot=annot.nki, do.mapping=TRUE,
  do.prediction.strength=FALSE, verbose=TRUE)
table(ssp2003.nki$subtype)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ clustering }
