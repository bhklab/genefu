\name{genius}
\alias{genius}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the Gene Expression progNostic Index Using Subtypes (GENIUS) as published by Haibe-Kains et al. 2010
}
\description{
This function computes the Gene Expression progNostic Index Using Subtypes (GENIUS) as published by Haibe-Kains et al. 2010. Subtype-specific risk scores are computed for each subtype signature separately and an overall risk score is computed by combining these scores with the posterior probability to belong to each of the breast cancer molecular subtypes.
}
\usage{
genius(data, annot, do.mapping = FALSE, mapping, do.scale = TRUE)
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
  \item{do.scale}{
\code{TRUE} if the ESR1, ERBB2 and AURKA (module) scores must be rescaled (see \code{\link[genefu]{rescale}}), \code{FALSE} otherwise.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{GENIUSM1 }{Risk score from the ER-/HER2- subtype signature in GENIUS model.}
  \item{GENIUSM2 }{Risk score from the HER2+ subtype signature in GENIUS model.}
  \item{GENIUSM3 }{Risk score from the ER+/HER2- subtype signature in GENIUS model.}
  \item{score }{Overall risk prediction as computed by the GENIUS model.}
}
\references{
Haibe-Kains B, Desmedt C, Rothe F, Sotiriou C and Bontempi G (2010) "A fuzzy gene expression-based computational approach improves breast cancer prognostication", \emph{Genome Biology}, \bold{11}(2):R18
}
\author{
Benjamin Haibe-Kains
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[genefu]{subtype.cluster.predict}},\code{\link[genefu]{sig.score}}
}
\examples{
## load NKI dataset
data(nki)
## compute GENIUS risk scores based on GENIUS model fitted on VDX dataset
genius.nki <- genius(data=data.nki, annot=annot.nki, do.mapping=TRUE)
str(genius.nki)
## the performance of GENIUS overall risk score predictions are not optimal
## since only part of the NKI dataset was used
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ prognosis }
