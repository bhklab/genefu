\name{claudinLowData}
\alias{claudinLowData}
\docType{data}
\title{
claudinLowData for use in the claudinLow classifier. Data generously provided by Aleix Prat. 
}
\description{
Training and Testing Data for use with the Claudin-Low Classifier
}
\usage{data(claudinLowData)}
%%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%%}
%\source{
%\url{http://jnci.oxfordjournals.org/cgi/content/full/98/4/262/DC1}
%}
\format{
\code{\link[genefu]{claudinLowData}} is a dataset containing one list with six sublists:
 \describe{
   \item{claudinLowData}{xd: Matrix of 807 features and 52 samples. classes: factor to split samples. nfeatures: number of features. nsamples: number of samples. fnames: names of features. snames: names of samples.}
}
}
\source{
  \url{http://www.biomedcentral.com/content/pdf/bcr2635.pdf}
}
\references{
Aleix Prat, Joel S Parker, Olga Karginova, Cheng Fan, Chad Livasy, Jason I Herschkowitz, Xiaping He, and Charles M. Perou (2010) "Phenotypic and molecular characterization of the claudin-low intrinsic subtype of breast cancer", \emph{Breast Cancer Research}, \bold{12}(5):R68
}
\seealso{
\code{\link[genefu]{claudinLow}}
}
\examples{
data(claudinLowData)
str(claudinLowData)
sapply(claudinLowData, head)
}
\keyword{data}
