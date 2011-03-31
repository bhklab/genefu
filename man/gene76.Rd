\name{gene76}
\alias{gene76}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the Relapse Score as published by Wang et al. 2005
}
\description{
This function computes signature scores and risk classifications from gene expression values following the algorithm used for the Relapse Score (GENE76) as published by Wang et al. 2005.
}
\usage{
gene76(data, er)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Matrix of gene expressions with samples in rows and probes in columns, dimnames being properly defined.
}
  \item{er}{
Vector containing the estrogen receptor (ER) status of breast cancer patients in the dataset.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
	\item{score }{Continuous signature scores}
	\item{risk }{Binary risk classification, \code{1} being high risk and \code{0} being low risk.}
}
\references{
Y. Wang and J. G. Klijn and Y. Zhang and A. M. Sieuwerts and M. P. Look and F. Yang and D. Talantov and M. Timmermans and M. E. Meijer-van Gelder and J. Yu and T. Jatkoe and E. M. Berns and D. Atkins and J. A. Foekens (2005) "Gene-Expression Profiles to Predict Distant Metastasis of Lymph-Node-Negative Primary Breast Cancer", \emph{Lancet}, \bold{365}(9460):671--679.
}
\author{
Benjamin Haibe-Kains
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[genefu]{ggi}}
}
\examples{
## load GENE76 signature
data(sig.gene76)
## load VDX dataset
data(vdx)
## compute relapse score
rs.vdx <- gene76(data=data.vdx, er=demo.vdx[ ,"er"])
table(rs.vdx$risk)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ prognosis }
%%\keyword{ breast cancer }
