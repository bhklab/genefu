\name{ihc4}
\alias{ihc4}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the IHC4 prognostic score as published by Paik et al. in 2004.
}
\description{
This function computes the prognostic score based on four measured IHC markers (ER, PGR, HER2, Ki-67), 
following the algorithm as published by Cuzick et al. 2011. The user has the option to either obtain just the 
shrinkage-adjusted IHC4 score (IHC4) or the overall score htat also combines the clinical score (IHC4+C)
}
\usage{
ihc4(ER, PGR, HER2, Ki67,age,size,grade,node,ana,scoreWithClinical=FALSE, na.rm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ER}{ER score between 0-10, calculated as (H-score/30)}
  \item{PGR}{Progesterone Receptor score between 0-10}
  \item{HER2}{Her2/neu status (0 or 1).}
  \item{Ki67}{Ki67 score based on percentage of positively staining malignant cells}
  \item{age}{patient age}
  \item{size}{tumor size in cm.}
  \item{grade}{Histological grade, i.e. low (1), intermediate (2) and high (3) grade.}
  \item{node}{Nodal status. }
  \item{ana}{treatment with anastrozole}
  \item{scoreWithClinical}{\code{TRUE} to get IHC4+C score, \code{FALSE} to get just the IHC4 score.}
  \item{na.rm}{\code{TRUE} if missing values should be removed, \code{FALSE} otherwise.}
}

%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}

\value{
Shrinkage-adjusted IHC4 score or the Overall Prognostic Score based on IHC4+C (IHC4+Clinical Score)
}

\references{
Jack Cuzick, Mitch Dowsett, Silvia Pineda, Christopher Wale, Janine Salter, Emma Quinn, Lila Zabaglo, Elizabeth Mallon, Andrew R. Green, Ian O. Ellis, Anthony Howell, Aman U. Buzdar, and John F. Forbes (2011) "Prognostic Value of a Combined Estrogen Receptor, Progesterone Receptor, Ki-67, and Human Epidermal Growth Factor Receptor 2 Immunohistochemical Score and Comparison with the Genomic Health Recurrence Score in Early Breast Cancer", \emph{Journal of Clinical Oncologoy}, \bold{29}(32):4273--4278.
}

\author{
Deena M.A. Gendoo
}
%%\note{
%%  ~~further notes~~
%%}

\examples{
## load NKI dataset
data(nkis)
## compute shrinkage-adjusted IHC4 score
count<-nrow(demo.nkis)
ihc4(ER=sample(x=1:10, size=count,replace=TRUE),PGR=sample(x=1:10, size=count,replace=TRUE),
HER2=sample(x=0:1,size=count,replace=TRUE),Ki67=sample(x=1:100, size=count,replace=TRUE),
scoreWithClinical=FALSE, na.rm=TRUE)

## compute IHC4+C score
ihc4(ER=sample(x=1:10, size=count,replace=TRUE),PGR=sample(x=1:10, size=count,replace=TRUE),
HER2=sample(x=0:1,size=count,replace=TRUE),Ki67=sample(x=1:100, size=count,replace=TRUE),
age=demo.nkis[,"age"],size=demo.nkis[ ,"size"],grade=demo.nkis[ ,"grade"],node=demo.nkis[ ,"node"],
ana=sample(x=0:1,size=count,replace=TRUE), scoreWithClinical=TRUE, na.rm=TRUE)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ prognosis }
%%\keyword{ breast cancer }

data(NKI)

