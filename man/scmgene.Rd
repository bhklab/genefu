\name{scmgene}
\alias{scmgene}
\docType{data}
\title{
Subtype Clustering Model using only ESR1, ERBB2 and AURKA genes for identification of breast cancer molecular subtypes
}
\description{
List of parameters defining the Subtype Clustering Model as published in Wirapati et al 2009 and Desmedt et al 2008 but using single genes instead of gene modules. 
}
\usage{data(scmgene)}
\format{
List of parameters for SCMGENE:
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
\url{http://clincancerres.aacrjournals.org/content/14/16/5158.abstract?ck=nck}
}
\references{
Desmedt C, Haibe-Kains B, Wirapati P, Buyse M, Larsimont D, Bontempi G, Delorenzi M, Piccart M, and Sotiriou C (2008) "Biological processes associated with breast cancer clinical outcome depend on the molecular subtypes", \emph{Clinical Cancer Research}, \bold{14}(16):5158--5165.
}
\examples{
data(scmgene)
str(scmgene, max.level=1)
}
\keyword{data}
