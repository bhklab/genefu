% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{scmod2.robust}
\alias{scmod2.robust}
\title{Subtype Clustering Model using ESR1, ERBB2 and AURKA modules for identification of breast cancer molecular subtypes (Desmedt et al 2008)}
\format{
List of parameters for SCMOD2:
\itemize{
\item parameters: List of parameters for the mixture of three Gaussians (ER-/HER2-, HER2+ and ER+/HER2-) that define the Subtype Clustering Model. The structure is the same than for an \code{\link[mclust:Mclust]{mclust::Mclust}} object.
\item cutoff.AURKA: Cutoff for AURKA module score in order to identify ER+/HER2- High Proliferation (aka Luminal B) tumors and ER+/HER2- Low Proliferation (aka Luminal A) tumors.
\item mod: ESR1, ERBB2 and AURKA modules.
}
}
\source{
\url{http://breast-cancer-research.com/content/10/4/R65k}
}
\usage{
data(scmod1.robust)
}
\description{
List of parameters defining the Subtype Clustering Model as published in Desmedt et al 2008.
}
\references{
Wirapati P, Sotiriou C, Kunkel S, Farmer P, Pradervand S, Haibe-Kains B, Desmedt C, Ignatiadis M, Sengstag T, Schutz F, Goldstein DR, Piccart MJ and Delorenzi M (2008) "Meta-analysis of Gene-Expression Profiles in Breast Cancer: Toward a Unified Understanding of Breast Cancer Sub-typing and Prognosis Signatures", Breast Cancer Research, 10(4):R65.
}
