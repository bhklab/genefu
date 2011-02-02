\name{vdx}
\alias{vdx}
\alias{data.vdx}
\alias{annot.vdx}
\alias{demo.vdx}
\docType{data}
\title{
Gene expression, annotations and clinical data from Wang et al. 2005 and Minn et al 2007
}
\description{
This dataset contains (part of) the gene xpression, annotations and clinical data as published in Wang et al. 2005 and Minn et al 2007.
}
\usage{data(vdx)}
\format{
\code{\link[genefu]{vdx}} is a dataset containing three matrices:
 \describe{
    \item{data.vdx}{Matrix containing gene expressions as measured by Affymetrix hgu133a technology (single-channel, oligonucleotides}
    \item{annot.vdx}{Matrix containing annotations of Agilent microarray platform}
    \item{demo.vdx}{Clinical information of the breast cancer patients whose tumors were hybridized}
}
}
\details{
This dataset represent only partially the one published by Wang et al. 2005 and Minn et al 2007. Indeed only part of the patients (150) and gene expressions (914) are contained in \code{\link[genefu]{data.vdx}}.
}
\source{
\url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2034}
\url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5327}
}
\references{
Y. Wang and J. G. Klijn and Y. Zhang and A. M. Sieuwerts and M. P. Look and F. Yang and D. Talantov and M. Timmermans and M. E. Meijer-van Gelder and J. Yu and T. Jatkoe and E. M. Berns and D. Atkins and J. A. Foekens (2005) "Gene-Expression Profiles to Predict Distant Metastasis of Lymph-Node-Negative Primary Breast Cancer", \emph{Lancet}, \bold{365}:671--679

Minn, Andy J. and Gupta, Gaorav P. and Padua, David and Bos, Paula and Nguyen, Don X. and Nuyten, Dimitry and Kreike, Bas and Zhang, Yi and Wang, Yixin and Ishwaran, Hemant and Foekens, John A. and van de Vijver, Marc and Massague, Joan (2007) "Lung metastasis genes couple breast tumor size and metastatic spread", \emph{Proceedings of the National Academy of Sciences}, \bold{104}(16):6740--6745
}
\examples{
data(vdx)
}
\keyword{data}
