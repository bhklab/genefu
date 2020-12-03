\name{molecular.subtyping}
\alias{molecular.subtyping}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to identify breast cancer molecular subtypes using the Subtype Clustering Model
}
\description{
This function identifies the breast cancer molecular subtypes using a Subtype Clustering Model fitted by \code{\link[genefu]{subtype.cluster}}.
}
\usage{
molecular.subtyping(sbt.model = c("scmgene", "scmod1", "scmod2", 
  "pam50", "ssp2006", "ssp2003", "intClust", "AIMS","claudinLow"), 
  data, annot, do.mapping = FALSE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sbt.model}{
Subtyping classification model, can be either "scmgene", "scmod1", "scmod2", "pam50", "ssp2006", "ssp2003", "intClust", "AIMS", or "claudinLow".
}
  \item{data}{
Matrix of gene expressions with samples in rows and probes in columns, dimnames being properly defined.
}
  \item{annot}{
Matrix of annotations with at least one column named "EntrezGene.ID" (for ssp, scm, AIMS, and claudinLow models) or "Gene.Symbol" (for the intClust model), dimnames being properly defined.
}
  \item{do.mapping}{
\code{TRUE} if the mapping through Entrez Gene ids must be performed (in case of ambiguities, the most variant probe is kept for each gene), \code{FALSE} otherwise.
}
  \item{verbose}{
\code{TRUE} if informative messages should be displayed, \code{FALSE} otherwise.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
  \item{subtype }{Subtypes identified by the subtyping classification model.}
  \item{subtype.proba }{Probabilities to belong to each subtype estimated by the subtyping classification model.}
  \item{subtype.crisp }{Crisp classes identified by the subtyping classification model.}
}
\references{
T. Sorlie and R. Tibshirani and J. Parker and T. Hastie and J. S. Marron and A. Nobel and S. Deng and H. Johnsen and R. Pesich and S. Geister and J. Demeter and C. Perou and P. E. Lonning and P. O. Brown and A. L. Borresen-Dale and D. Botstein (2003) "Repeated Observation of Breast Tumor Subtypes in Independent Gene Expression Data Sets", \emph{Proceedings of the National Academy of Sciences}, \bold{1}(14):8418--8423

Hu, Zhiyuan and Fan, Cheng and Oh, Daniel and Marron, JS and He, Xiaping and Qaqish, Bahjat and Livasy, Chad and Carey, Lisa and Reynolds, Evangeline and Dressler, Lynn and Nobel, Andrew and Parker, Joel and Ewend, Matthew and Sawyer, Lynda and Wu, Junyuan and Liu, Yudong and Nanda, Rita and Tretiakova, Maria and Orrico, Alejandra and Dreher, Donna and Palazzo, Juan and Perreard, Laurent and Nelson, Edward and Mone, Mary and Hansen, Heidi and Mullins, Michael and Quackenbush, John and Ellis, Matthew and Olopade, Olufunmilayo and Bernard, Philip and Perou, Charles (2006) "The molecular portraits of breast tumors are conserved across microarray platforms", \emph{BMC Genomics}, \bold{7}(96)

Parker, Joel S. and Mullins, Michael and Cheang, Maggie C.U. and Leung, Samuel and Voduc, David and Vickery, Tammi and Davies, Sherri and Fauron, Christiane and He, Xiaping and Hu, Zhiyuan and Quackenbush, John F. and Stijleman, Inge J. and Palazzo, Juan and Marron, J.S. and Nobel, Andrew B. and Mardis, Elaine and Nielsen, Torsten O. and Ellis, Matthew J. and Perou, Charles M. and Bernard, Philip S. (2009) "Supervised Risk Predictor of Breast Cancer Based on Intrinsic Subtypes", \emph{Journal of Clinical Oncology}, \bold{27}(8):1160--1167

Desmedt C, Haibe-Kains B, Wirapati P, Buyse M, Larsimont D, Bontempi G, Delorenzi M, Piccart M, and Sotiriou C (2008) "Biological processes associated with breast cancer clinical outcome depend on the molecular subtypes", \emph{Clinical Cancer Research}, \bold{14}(16):5158--5165.

Wirapati P, Sotiriou C, Kunkel S, Farmer P, Pradervand S, Haibe-Kains B, Desmedt C, Ignatiadis M, Sengstag T, Schutz F, Goldstein DR, Piccart MJ and Delorenzi M (2008) "Meta-analysis of Gene-Expression Profiles in Breast Cancer: Toward a Unified Understanding of Breast Cancer Sub-typing and Prognosis Signatures", \emph{Breast Cancer Research}, \bold{10}(4):R65.

Haibe-Kains B, Desmedt C, Loi S, Culhane AC, Bontempi G, Quackenbush J, Sotiriou C. (2012) "A three-gene model to robustly identify breast cancer molecular subtypes.", \emph{J Natl Cancer Inst.}, \bold{104}(4):311--325.

Curtis C, Shah SP, Chin SF, Turashvili G, Rueda OM, Dunning MJ, Speed D, Lynch AG, Samarajiwa S, Yuan Y, Graf S, Ha G, Haffari G, Bashashati A, Russell R, McKinney S; METABRIC Group, Langerod A, Green A, Provenzano E, Wishart G, Pinder S, Watson P, Markowetz F, Murphy L, Ellis I, Purushotham A, Borresen-Dale AL, Brenton JD, Tavare S, Caldas C, Aparicio S. (2012) "The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups.", \emph{Nature}, \bold{486}(7403):346--352.

Paquet ER, Hallett MT. (2015) "Absolute assignment of breast cancer intrinsic molecular subtype.", \emph{J Natl Cancer Inst.}, \bold{107}(1):357.

Aleix Prat, Joel S Parker, Olga Karginova, Cheng Fan, Chad Livasy, Jason I Herschkowitz, Xiaping He, and Charles M. Perou (2010) "Phenotypic and molecular characterization of the claudin-low intrinsic subtype of breast cancer", \emph{Breast Cancer Research}, \bold{12}(5):R68

}
\author{
Benjamin Haibe-Kains
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[genefu]{subtype.cluster.predict}}, \code{\link[genefu]{intrinsic.cluster.predict}}
}
\examples{
##### without mapping (affy hgu133a or plus2 only)
## load VDX data
data(vdxs)

## Subtype Clustering Model fitted on EXPO and applied on VDX
sbt.vdx.SCMGENE <- molecular.subtyping(sbt.model="scmgene", 
  data=data.vdxs, annot=annot.vdxs, do.mapping=FALSE)
table(sbt.vdx.SCMGENE$subtype)

## Using the AIMS molecular subtyping algorithm
sbt.vdxs.AIMS <- molecular.subtyping(sbt.model="AIMS", data=data.vdxs,
  annot=annot.vdxs, do.mapping=FALSE)
table(sbt.vdxs.AIMS$subtype)

## Using the IntClust molecular subtyping algorithm
colnames(annot.vdxs)[3]<-"Gene.Symbol"
sbt.vdxs.intClust <- molecular.subtyping(sbt.model="intClust", data=data.vdxs,
  annot=annot.vdxs, do.mapping=FALSE)
table(sbt.vdxs.intClust$subtype)

##### with mapping
## load NKI data
data(nkis)

## Subtype Clustering Model fitted on EXPO and applied on NKI
sbt.nkis <- molecular.subtyping(sbt.model="scmgene", data=data.nkis, 
  annot=annot.nkis, do.mapping=TRUE)
table(sbt.nkis$subtype)

##### with mapping
## load vdxs data
data(vdxs)

## Claudin-Low classification of 150 VDXS samples
sbt.vdxs.CL <- molecular.subtyping(sbt.model="claudinLow", data=data.vdxs, 
  annot=annot.vdxs, do.mapping=TRUE)
table(sbt.vdxs.CL$subtype)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ clustering }
