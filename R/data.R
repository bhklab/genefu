#' @name claudinLowData
#'
#' @title claudinLowData for use in the claudinLow classifier. Data generously provided by Aleix Prat.
#'
#' @description Training and Testing Data for use with the Claudin-Low Classifier
#'
#' @usage data(claudinLowData)
#'
#' @source http://jnci.oxfordjournals.org/cgi/content/full/98/4/262/DC1
#'
#' @format
#'  - xd: Matrix of 807 features and 52 samples
#'  - classes: factor to split samples
#'  - nfeatures: number of features
#'  - nsamples: number of samples
#'  - fnames: names of features
#'  - snames: names of samples
#'
#' @references Aleix Prat, Joel S Parker, Olga Karginova, Cheng Fan, Chad Livasy, Jason I Herschkowitz, Xiaping He, and Charles M. Perou (2010) "Phenotypic and molecular characterization of the claudin-low intrinsic subtype of breast cancer", Breast Cancer Research, 12(5):R68
#'
#' @seealso [genefu::claudinLow]
#'
#' @md
#' @docType data
#' @keywords data
NULL

#' @name expos
#'
#' @aliases data.expos annot.expos demo.expos
#'
#' @title Gene expression, annotations and clinical data from the International Genomics Consortium
#'
#' @description This dataset contains (part of) the gene expression, annotations and clinical data from the expO dataset collected by the International Genomics Consortium ([](http://www.intgen.org/expo/)).
#'
#' @format expos is a dataset containing three matrices
#'   - data.expos: Matrix containing gene expressions as measured by Affymetrix hgu133plus2 technology (single-channel, oligonucleotides)
#'   - annot.expos: Matrix containing annotations of ffymetrix hgu133plus2 microarray platform
#'   - demo.expos: Clinical information of the breast cancer patients whose tumors were hybridized
#'
#' @source http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2109
#'
#' @references
#' International Genomics Consortium, http://www.intgen.org/research-services/biobanking-experience/expo/
#' McCall MN, Bolstad BM, Irizarry RA. (2010) "Frozen robust multiarray analysis (fRMA)", Biostatistics, 11(2):242-253.
#'
#' @usage data(expos)
#'
#' @md
#' @docType data
#' @keywords data
NULL

#' @name mod1
#' @md
#' @docType data
#' @title Gene modules published in Desmedt et al. 2008
#' @description List of seven gene modules published in Desmedt et a. 2008, i.e. ESR1 (estrogen receptor pathway), ERBB2 (her2/neu receptor pathway), AURKA (proliferation), STAT1 (immune response), PLAU (tumor invasion), VEGF (angogenesis) and CASP3 (apoptosis).
#' @usage data(mod1)
#' @details mod1 is a list of seven gene signatures, i.e. matrices with 3 columns containing the annotations and information related to the signatures themselves.
#' @references Desmedt C, Haibe-Kains B, Wirapati P, Buyse M, Larsimont D, Bontempi G, Delorenzi M, Piccart M, and Sotiriou C (2008) "Biological processes associated with breast cancer clinical outcome depend on the molecular subtypes", Clinical Cancer Research, 14(16):5158--5165.
#' @keywords data
NULL

#' @name mod2
#' @md
#' @docType data
#' @title Gene modules published in Wirapati et al. 2008
#' @description List of seven gene modules published in Wirapati et a. 2008, i.e. ESR1 (estrogen receptor pathway), ERBB2 (her2/neu receptor pathway) and AURKA (proliferation).
#' @usage data(mod2)
#' @details mod2 is a list of three gene signatures, i.e. matrices with 3 columns containing the annotations and information related to the signatures themselves.
#' @source [](http://breast-cancer-research.com/content/10/4/R65)
#' @references Wirapati P, Sotiriou C, Kunkel S, Farmer P, Pradervand S, Haibe-Kains B, Desmedt C, Ignatiadis M, Sengstag T, Schutz F, Goldstein DR, Piccart MJ and Delorenzi M (2008) "Meta-analysis of Gene-Expression Profiles in Breast Cancer: Toward a Unified Understanding of Breast Cancer Sub-typing and Prognosis Signatures", Breast Cancer Research, 10(4):R65.
#' @keywords data
NULL

#' @name modelOvcAngiogenic
#' @docType data
#' @title Model used to classify ovarian tumors into Angiogenic and NonAngiogenic subtypes.
#' @description Object containing the set of parameters for the mixture of Gaussians used as a model to classify ovarian tumors into Angiogenic and NonAngiogenic subtypes.
#' @usage data(modelOvcAngiogenic)
#' @source [](http://jnci.oxfordjournals.org/cgi/content/full/98/4/262/DC1)
#' @references Bentink S, Haibe-Kains B, Risch T, Fan J-B, Hirsch MS, Holton K, Rubio R, April C, Chen J, Wickham-Garcia E, Liu J, Culhane AC, Drapkin R, Quackenbush JF, Matulonis UA (2012) "Angiogenic mRNA and microRNA Gene Expression Signature Predicts a Novel Subtype of Serous Ovarian Cancer", PloS one, 7(2):e30269
#' @keywords data
NULL

#' @name nkis
#' @md
#' @aliases data.nkis annot.nkis demo.nkis
#' @docType data
#' @title Gene expression, annotations and clinical data from van de Vijver et al. 2002
#' @description This dataset contains (part of) the gene expression, annotations and clinical data as published in van de Vijver et al. 2002.
#' @usage data(nkis)
#' @details This dataset represent only partially the one published by van  de Vijver et al. in 2008. Indeed, only part of the patients (150) and gene expressions (922) in [`data.nkis`].
#' @format nkis is a dataset containing three matrices:
#'   - data.nkis: Matrix containing gene expressions as measured by Agilent technology (dual-channel, oligonucleotides)
#'   - annot.nkis: Matrix containing annotations of Agilent microarray platform
#'   - demon.nkis: Clinical information of the breast cancer patients whose tumors were hybridized
#' @source http://www.nature.com/nature/journal/v415/n6871/full/415530a.html
#' @references M. J. van de Vijver and Y. D. He and L. van't Veer and H. Dai and A. M. Hart and D. W. Voskuil and G. J. Schreiber and J. L. Peterse and C. Roberts and M. J. Marton and M. Parrish and D. Atsma and A. Witteveen and A. Glas and L. Delahaye and T. van der Velde and H. Bartelink and S. Rodenhuis and E. T. Rutgers and S. H. Friend and R. Bernards (2002) "A Gene Expression Signature as a Predictor of Survival in Breast Cancer", New England Journal of Medicine, 347(25):1999--2009
#' @keywords data
NULL

#' @name pam50
#' @md
#' @aliases pam50.scale pam50.robust
#' @docType data
#' @title PAM50 classifier for identification of breast cancer molecular subtypes (Parker et al 2009)
#' @description List of parameters defining the PAM50 classifier for identification of breast cancer molecular subtypes (Parker et al 2009).
#' @usage
#' data(pam50)
#' data(pam50.scale)
#' data(pam50.robust)
#' @format List of parameters for PAM50:
#'  - centroids: Gene expression centroids for each subtype.
#'  - centroids.map: Mapping for centroids.
#'  - method.cor: Method of correlation used to compute distance to the centroids.
#'  - method.centroids: Method used to compute the centroids.
#'  - std: Method of standardization for gene expressions ("none", "scale" or "robust")
#'  - mins: Minimum number of samples within each cluster allowed during the fitting of the model.
#' @details Three versions of the model are provided, each of ones differs by the gene expressions standardization method since it has an important impact on the subtype classification:
#'   - pam50: Use of the official centroids without scaling of the gene expressions.
#'   - pam50.scale: Use of the official centroids with traditional scaling of the gene expressions (see [`base::scale`])
#'   - pam50.robust: Use of the official centroids with robust scaling of the gene expressions (see [`genefu::rescale`])
#' The model `pam50.robust`` has been shown to reach the best concordance with the traditional clinical parameters (ER IHC, HER2 IHC/FISH and histological grade). However the use of this model is recommended only when the dataset is representative of a global population of breast cancer patients (no sampling bias, the 5 subtypes should be present).
#' @source http://jco.ascopubs.org/cgi/content/short/JCO.2008.18.1370v1
#' @references Parker, Joel S. and Mullins, Michael and Cheang, Maggie C.U. and Leung, Samuel and Voduc, David and Vickery, Tammi and Davies, Sherri and Fauron, Christiane and He, Xiaping and Hu, Zhiyuan and Quackenbush, John F. and Stijleman, Inge J. and Palazzo, Juan and Marron, J.S. and Nobel, Andrew B. and Mardis, Elaine and Nielsen, Torsten O. and Ellis, Matthew J. and Perou, Charles M. and Bernard, Philip S. (2009) "Supervised Risk Predictor of Breast Cancer Based on Intrinsic Subtypes", Journal of Clinical Oncology, 27(8):1160--1167
#' @keywords data
NULL

#' @name scmgene.robust
#' @docType data
#' @title Subtype Clustering Model using only ESR1, ERBB2 and AURKA genes for identification of breast cancer molecular subtypes
#' @description List of parameters defining the Subtype Clustering Model as published in Wirapati et al 2009 and Desmedt et al 2008 but using single genes instead of gene modules.
#' @usage data(scmgene.robust)
#' @format List of parameters for SCMGENE:
#'   - parameters: List of parameters for the mixture of three Gaussians (ER-/HER2-, HER2+ and ER+/HER2-) that define the Subtype Clustering Model. The structure is the same than for an [`mclust::Mclust`] object.
#'   - cutoff.AURKA: Cutoff for AURKA module score in order to identify ER+/HER2- High Proliferation (aka Luminal B) tumors and ER+/HER2- Low Proliferation (aka Luminal A) tumors.
#'   - mod: ESR1, ERBB2 and AURKA modules.
#' @source http://clincancerres.aacrjournals.org/content/14/16/5158.abstract?ck=nck
#' @references Desmedt C, Haibe-Kains B, Wirapati P, Buyse M, Larsimont D, Bontempi G, Delorenzi M, Piccart M, and Sotiriou C (2008) "Biological processes associated with breast cancer clinical outcome depend on the molecular subtypes", Clinical Cancer Research, 14(16):5158--5165.
#' @keywords data
#' @md
NULL

#' @name scmod1.robust
#' @docType data
#' @title Subtype Clustering Model using ESR1, ERBB2 and AURKA modules for identification of breast cancer molecular subtypes (Desmedt et al 2008)
#' @description List of parameters defining the Subtype Clustering Model as published in Desmedt et al 2008.
#' @usage data(scmod1.robust)
#' @format List of parameters for SCMOD1:
#'   - parameters: List of parameters for the mixture of three Gaussians (ER-/HER2-, HER2+ and ER+/HER2-) that define the Subtype Clustering Model. The structure is the same than for an [`mclust::Mclust`] object.
#'   - cutoff.AURKA: Cutoff for AURKA module score in order to identify ER+/HER2- High Proliferation (aka Luminal B) tumors and ER+/HER2- Low Proliferation (aka Luminal A) tumors.
#'   - mod: ESR1, ERBB2 and AURKA modules.
#' @source [](http://clincancerres.aacrjournals.org/content/14/16/5158.abstract?ck=nck)
#' @references Desmedt C, Haibe-Kains B, Wirapati P, Buyse M, Larsimont D, Bontempi G, Delorenzi M, Piccart M, and Sotiriou C (2008) "Biological processes associated with breast cancer clinical outcome depend on the molecular subtypes", Clinical Cancer Research}, 14(16):5158--5165.
#' @keywords data
#' @md
NULL

#' @name scmod2.robust
#' @docType data
#' @title Subtype Clustering Model using ESR1, ERBB2 and AURKA modules for identification of breast cancer molecular subtypes (Desmedt et al 2008)
#' @description List of parameters defining the Subtype Clustering Model as published in Desmedt et al 2008.
#' @usage data(scmod1.robust)
#' @format List of parameters for SCMOD2:
#'   - parameters: List of parameters for the mixture of three Gaussians (ER-/HER2-, HER2+ and ER+/HER2-) that define the Subtype Clustering Model. The structure is the same than for an [`mclust::Mclust`] object.
#'   - cutoff.AURKA: Cutoff for AURKA module score in order to identify ER+/HER2- High Proliferation (aka Luminal B) tumors and ER+/HER2- Low Proliferation (aka Luminal A) tumors.
#'   - mod: ESR1, ERBB2 and AURKA modules.
#' @source [](http://breast-cancer-research.com/content/10/4/R65k)
#' @references Wirapati P, Sotiriou C, Kunkel S, Farmer P, Pradervand S, Haibe-Kains B, Desmedt C, Ignatiadis M, Sengstag T, Schutz F, Goldstein DR, Piccart MJ and Delorenzi M (2008) "Meta-analysis of Gene-Expression Profiles in Breast Cancer: Toward a Unified Understanding of Breast Cancer Sub-typing and Prognosis Signatures", Breast Cancer Research, 10(4):R65.
#' @md
NULL

#' @name sig.endoPredict
#' @docType data
#' @title Signature used to compute the endoPredict signature as published by Filipits et al 2011
#' @description List of 11 genes included in the endoPredict signature. The EntrezGene.ID allows for mapping and the mapping to affy probes is already provided.
#' @usage data(sig.endoPredict)
#' @format sig.endoPredict is a matrix with 5 columns containing the annotations and information related to the signature itself (including a mapping to Affymetrix HGU platform).
#' @reference Filipits, M., Rudas, M., Jakesz, R., Dubsky, P., Fitzal, F., Singer, C. F., et al. (2011). "A new molecular predictor of distant recurrence in ER-positive, HER2-negative breast cancer adds independent information to conventional clinical risk factors." \emph{Clinical Cancer Research}, \bold{17}(18):6012--6020.
#' @keywords data
#' @md
NULL

#' @name sig.gene70
#' @docType data
#' @title Signature used to compute the 70 genes prognosis profile (GENE70) as published by van't Veer et al. 2002
#' @description List of 70 agilent probe ids representing 56 unique genes included in the GENE70 signature. The EntrezGene.ID allows for mapping and the "average.good.prognosis.profile" values allows for signature computation.
#' @usage data(sig.gene70)
#' @format sig.gene70 is a matrix with 9 columns containing the annotations and information related to the signature itself.
#' @source [](http://www.nature.com/nature/journal/v415/n6871/full/415530a.html)
#' @references L. J. van't Veer and H. Dai and M. J. van de Vijver and Y. D. He and A. A. Hart and M. Mao and H. L. Peterse and K. van der Kooy and M. J. Marton and A. T. Witteveen and G. J. Schreiber and R. M. Kerkhiven and C. Roberts and P. S. Linsley and R. Bernards and S. H. Friend (2002) "Gene Expression Profiling Predicts Clinical Outcome of Breast Cancer", Nature, 415:530--536.
#' @keywords data
#' @md
NULL

#














