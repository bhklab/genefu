\name{scmod2}
\alias{scmod2}
\docType{data}
\title{
Subtype Clustering Model using ESR1, ERBB2 and AURKA modules for identification of breast cancer molecular subtypes (Wirapati et al 2008)
}
\description{
List of parameters defining the Subtype Clustering Model as published in Wirapati et al 2008. 
}
\usage{data(scmod2)}
\format{
List of parameters for SCMOD2:
\describe{
\item{\code{parameters} }{List of parameters for the mixture of three Gaussians (ER-/HER2-, HER2+ and ER+/HER2-) that define the Subtype Clustering Model. The structure is the same than for an \code{\link[mclust]{Mclust}} object.}
\item{\code{cutoff.AURKA} }{Cutoff for AURKA module score in order to identify ER+/HER2- High Proliferation (aka Luminal B) tumors and ER+/HER2- Low Proliferation (aka Luminal A) tumors. }
\item{\code{mod} }{ESR1, ERBB2 and AURKA modules.}
}
}
%%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%%}
\source{
\url{http://breast-cancer-research.com/content/10/4/R65}
}
\references{
Wirapati P, Sotiriou C, Kunkel S, Farmer P, Pradervand S, Haibe-Kains B, Desmedt C, Ignatiadis M, Sengstag T, Schutz F, Goldstein DR, Piccart MJ and Delorenzi M (2008) "Meta-analysis of Gene-Expression Profiles in Breast Cancer: Toward a Unified Understanding of Breast Cancer Sub-typing and Prognosis Signatures", \emph{Breast Cancer Research}, \bold{10}(4):R65.
}
\examples{
data(scmod2)
str(scmod2, max.level=1)
}
\keyword{data}
