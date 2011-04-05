\name{st.gallen}
\alias{st.gallen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the St Gallen consensus criterion for prognostication
}
\description{
This function computes the updated St Gallen consensus criterions as published by Goldhirsh et al 2003.
}
\usage{
st.gallen(size, grade, node, her2.neu, age, vascular.inv, na.rm = FALSE)
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
Nodal status (0 or 1 for no lymph node invasion a,d at least 1 invaded lymph ode respectively).
}
  \item{her2.neu}{
Her2/neu status (0 or 1).
}
  \item{age}{
Age at diagnosis (in years).
}
  \item{vascular.inv}{
Peritumoral vascular invasion (0 or 1).
}
  \item{na.rm}{
\code{TRUE} if missing values should be removed, \code{FALSE} otherwise.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
Vector of risk predictions: "Good", "Intermediate", and "Poor".
}
\references{
Goldhirsh A, Wood WC, Gelber RD, Coates AS, Thurlimann B, and Senn HJ (2003) "Meeting highlights: Updated international expert consensus on the primary therapy of early breast cancer", \emph{Journal of Clinical Oncology}, \bold{21}(17):3357--3365.
}
\author{
Benjamin Haibe-Kains
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[genefu]{npi}}
}
\examples{
## load NKI dataset
data(NKI)
## compute St Gallen predictions
st.gallen(size=demo.nki[ ,"size"], grade=demo.nki[ ,"grade"],
  node=demo.nki[ ,"node"], her2.neu=sample(x=0:1, size=nrow(demo.nki),
  replace=TRUE), age=demo.nki[ ,"age"], vascular.inv=sample(x=0:1,
  size=nrow(demo.nki), replace=TRUE), na.rm=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ prognosis }
