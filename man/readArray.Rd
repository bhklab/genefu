\name{readArray}
\alias{readArray}
\title{Overlap two datasets}
\description{
 Formatting function to read arrays and format for use in the claudinLow classifier
}
\usage{
readarray<-function(dataFile,designFile=NA,hr=1,impute=T,method="mean")
}
\arguments{
  \item{dataFile}{file with matrix to be read}
  \item{hr}{Header rows as Present (2) or Absent (1)}
  \item{designFile}{Design of file}
  \item{impute}{whether data will be imputed or not.}
}
\references{
  \code{citation("claudinLow")}
}
\author{}
\seealso{ \code{\link[claudinLow]{claudinLow}} }


