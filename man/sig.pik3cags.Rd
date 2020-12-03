\name{sig.pik3cags}
\alias{sig.pik3cags}
\docType{data}
\title{
Gene expression Grade Index (GGI) as published in Sotiriou et al. 2006
}
\description{
List of 278 affymetrix hgu133a probesets representing 236 unique genes included in the PIK3CA-GS signature. The "EntrezGene.ID" column allows for mapping and "coefficient" refers to to the direction of association with PIK3CA mutation.
}
\usage{data(sig.pik3cags)}
\format{
  \code{\link[genefu]{sig.pik3cags}} is a matrix with 3 columns containing the annotations and information related to the signature itself.
}
%%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%%}
\source{
\url{http://www.pnas.org/content/107/22/10208/suppl/DCSupplemental}
}
\references{
Loi S, Haibe-Kains B, Majjaj S, Lallemand F, Durbecq V, Larsimont D, Gonzalez-Angulo AM, Pusztai L, Symmans FW, Bardelli A, Ellis P, Tutt AN, Gillett CE, Hennessy BT., Mills GB, Phillips WA, Piccart MJ, Speed TP, McArthur GA, Sotiriou C (2010) "PIK3CA mutations associated with gene signature of low mTORC1 signaling and better outcomes in estrogen receptor-positive breast cancer", \emph{Proceedings of the National Academy of Sciences}, \bold{107}(22):10208--10213
}
\examples{
data(sig.pik3cags)
head(sig.pik3cags)
}
\keyword{data}
