\name{npi}
\alias{npi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the Nottingham Prognostic Index
}
\description{
This function computes the Nottingham Prognostic Index (NPI) as published in Galeat et al, 1992. NPI is a clinical index shown to be highly prognostic in breast cancer.
}
\usage{
npi(size, grade, node, na.rm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{size}{
tumor size in cm.
}
  \item{grade}{
Histological grade, i.e. low (1), intermediate (2) and high (3) grade.
}
  \item{node}{
Nodal status. If only binary nodal status (0/1) is available, map \code{0} to \code{1} and \code{1} to \code{3}.
}
  \item{na.rm}{
\code{TRUE} if missing values should be removed, \code{FALSE} otherwise.
}
}
\details{
The risk prediction is either \code{Good} if score < 3.4, \code{Intermediate} if 3.4 <= score <- 5.4, or \code{Poor} if score > 5.4.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{score }{Continuous signature scores}
\item{risk }{Binary risk classification, \code{1} being high risk and \code{0} being low risk.}
}
\references{
Galea MH, Blamey RW, Elston CE, and Ellis IO (1992) "The nottingham prognostic index in primary breast cancer", \emph{Breast Cancer Reasearch and Treatment}, \bold{22}(3):207--219.
}
\author{
Benjamin Haibe-Kains
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[genefu]{st.gallen}}
}
\examples{
## load NKI dataset
data(nki)
## compute NPI score and risk classification
npi(size=demo.nki[ ,"size"], grade=demo.nki[ ,"grade"],
  node=ifelse(demo.nki[ ,"node"] == 0, 1, 3), na.rm=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ prognosis }
