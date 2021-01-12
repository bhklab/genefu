#' @name claudinLowData
#'
#' @title claudinLowData for use in the claudinLow classifier. Data generously provided by Aleix Prat.
#'
#' @description Training and Testing Data for use with the Claudin-Low Classifier
#'
#' @usage data(claudinLowData)
#'
#' @source [http://jnci.oxfordjournals.org/cgi/content/full/98/4/262/DC1](http://jnci.oxfordjournals.org/cgi/content/full/98/4/262/DC1)
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
#' @seealso [genefu::claudinLow()]
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
#' @source [http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2109](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2109)
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
#' @source [http://breast-cancer-research.com/content/10/4/R65](http://breast-cancer-research.com/content/10/4/R65)
#' @references Wirapati P, Sotiriou C, Kunkel S, Farmer P, Pradervand S, Haibe-Kains B, Desmedt C, Ignatiadis M, Sengstag T, Schutz F, Goldstein DR, Piccart MJ and Delorenzi M (2008) "Meta-analysis of Gene-Expression Profiles in Breast Cancer: Toward a Unified Understanding of Breast Cancer Sub-typing and Prognosis Signatures", Breast Cancer Research, 10(4):R65.
#' @keywords data
NULL

#' @name modelOvcAngiogenic
#' @docType data
#' @title Model used to classify ovarian tumors into Angiogenic and NonAngiogenic subtypes.
#' @description Object containing the set of parameters for the mixture of Gaussians used as a model to classify ovarian tumors into Angiogenic and NonAngiogenic subtypes.
#' @usage data(modelOvcAngiogenic)
#' @source [http://jnci.oxfordjournals.org/cgi/content/full/98/4/262/DC1](http://jnci.oxfordjournals.org/cgi/content/full/98/4/262/DC1)
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
#' @source [http://www.nature.com/nature/journal/v415/n6871/full/415530a.html](http://www.nature.com/nature/journal/v415/n6871/full/415530a.html)
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
#'   - pam50.scale: Use of the official centroids with traditional scaling of the gene expressions (see [`base::scale()`])
#'   - pam50.robust: Use of the official centroids with robust scaling of the gene expressions (see [`genefu::rescale()`])
#' The model `pam50.robust`` has been shown to reach the best concordance with the traditional clinical parameters (ER IHC, HER2 IHC/FISH and histological grade). However the use of this model is recommended only when the dataset is representative of a global population of breast cancer patients (no sampling bias, the 5 subtypes should be present).
#' @source [http://jco.ascopubs.org/cgi/content/short/JCO.2008.18.1370v1](http://jco.ascopubs.org/cgi/content/short/JCO.2008.18.1370v1)
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
#' @source [http://clincancerres.aacrjournals.org/content/14/16/5158.abstract?ck=nck](http://clincancerres.aacrjournals.org/content/14/16/5158.abstract?ck=nck)
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
#'   - parameters: List of parameters for the mixture of three Gaussians (ER-/HER2-, HER2+ and ER+/HER2-) that define the Subtype Clustering Model. The structure is the same than for an [`mclust::Mclust()`] object.
#'   - cutoff.AURKA: Cutoff for AURKA module score in order to identify ER+/HER2- High Proliferation (aka Luminal B) tumors and ER+/HER2- Low Proliferation (aka Luminal A) tumors.
#'   - mod: ESR1, ERBB2 and AURKA modules.
#' @source [http://clincancerres.aacrjournals.org/content/14/16/5158.abstract?ck=nck](http://clincancerres.aacrjournals.org/content/14/16/5158.abstract?ck=nck)
#' @references Desmedt C, Haibe-Kains B, Wirapati P, Buyse M, Larsimont D, Bontempi G, Delorenzi M, Piccart M, and Sotiriou C (2008) "Biological processes associated with breast cancer clinical outcome depend on the molecular subtypes", _Clinical Cancer Research_, *14*(16):5158--5165.
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
#' @source [http://breast-cancer-research.com/content/10/4/R65k](http://breast-cancer-research.com/content/10/4/R65k)
#' @references Wirapati P, Sotiriou C, Kunkel S, Farmer P, Pradervand S, Haibe-Kains B, Desmedt C, Ignatiadis M, Sengstag T, Schutz F, Goldstein DR, Piccart MJ and Delorenzi M (2008) "Meta-analysis of Gene-Expression Profiles in Breast Cancer: Toward a Unified Understanding of Breast Cancer Sub-typing and Prognosis Signatures", Breast Cancer Research, 10(4):R65.
#' @md
NULL

#' @name sig.endoPredict
#' @docType data
#' @title Signature used to compute the endoPredict signature as published by Filipits et al 2011
#' @description List of 11 genes included in the endoPredict signature. The EntrezGene.ID allows for mapping and the mapping to affy probes is already provided.
#' @usage data(sig.endoPredict)
#' @format `sig.endoPredict` is a matrix with 5 columns containing the annotations and information related to the signature itself (including a mapping to Affymetrix HGU platform).
#' @references Filipits, M., Rudas, M., Jakesz, R., Dubsky, P., Fitzal, F., Singer, C. F., et al. (2011). "A new molecular predictor of distant recurrence in ER-positive, HER2-negative breast cancer adds independent information to conventional clinical risk factors." \emph{Clinical Cancer Research}, \bold{17}(18):6012--6020.
#' @keywords data
#' @md
NULL

#' @name sig.gene70
#' @docType data
#' @title Signature used to compute the 70 genes prognosis profile (GENE70) as published by van't Veer et al. 2002
#' @description List of 70 agilent probe ids representing 56 unique genes included in the GENE70 signature. The EntrezGene.ID allows for mapping and the "average.good.prognosis.profile" values allows for signature computation.
#' @usage data(sig.gene70)
#' @format sig.gene70 is a matrix with 9 columns containing the annotations and information related to the signature itself.
#' @source [http://www.nature.com/nature/journal/v415/n6871/full/415530a.html](http://www.nature.com/nature/journal/v415/n6871/full/415530a.html)
#' @references L. J. van't Veer and H. Dai and M. J. van de Vijver and Y. D. He and A. A. Hart and M. Mao and H. L. Peterse and K. van der Kooy and M. J. Marton and A. T. Witteveen and G. J. Schreiber and R. M. Kerkhiven and C. Roberts and P. S. Linsley and R. Bernards and S. H. Friend (2002) "Gene Expression Profiling Predicts Clinical Outcome of Breast Cancer", Nature, 415:530--536.
#' @keywords data
#' @md
NULL

#' @name sig.gene76
#' @docType data
#' @title Signature used to compute the Relapse Score (GENE76) as published in Wang et al. 2005
#' @description List of 76 affymetrix hgu133a probesets representing 60 unique genes included in the GENE76 signature. The EntrezGene.ID allows for mapping and the coefficient allows for signature computation.
#' @usage data(sig.gene76)
#' @format `sig.gene70` is a matrix with 10 columns containing the annotations and information related to the signature itself.
#' @source [http://www.thelancet.com/journals/lancet/article/PIIS0140-6736(05)17947-1/abstract](http://www.thelancet.com/journals/lancet/article/PIIS0140-6736(05)17947-1/abstract)
#' @references Y. Wang and J. G. Klijn and Y. Zhang and A. M. Sieuwerts and M. P. Look and F. Yang and D. Talantov and M. Timmermans and M. E. Meijer-van Gelder and J. Yu and T. Jatkoe and E. M. Berns and D. Atkins and J. A. Foekens (2005) "Gene-Expression Profiles to Predict Distant Metastasis of Lymph-Node-Negative Primary Breast Cancer", Lancet, 365(9460):671--679.
#' @keywords data
#' @md
NULL

#' @name sig.genius
#' @docType data
#' @title Gene Expression progNostic Index Using Subtypes (GENIUS) as published by Haibe-Kains et al. 2010.
#' @description List of three gene signatures which compose the Gene Expression progNostic Index Using Subtypes (GENIUS) as published by Haibe-Kains et al. 2009. GENIUSM1, GENIUSM2 and GENIUSM3  are the ER-/HER2-, HER2+ and ER+/HER2- subtype signatures respectively.
#' @format `sig.genius` is a list a three subtype signatures.
#' @references Haibe-Kains B, Desmedt C, Rothe F, Sotiriou C and Bontempi G (2010) "A fuzzy gene expression-based computational approach improves breast cancer prognostication", Genome Biology, 11(2):R18
#' @keywords data
#' @md
NULL

#' @name sig.ggi
#' @docType data
#' @title Gene expression Grade Index (GGI) as published in Sotiriou et al. 2006
#' @description List of 128 affymetrix hgu133a probesets representing 97 unique genes included in the GGI signature. The "EntrezGene.ID" column allows for mapping and "grade" defines the up-regulation of the expressions either in histological grade 1 or 3.
#' @usage data(sig.ggi)
#' @format sig.ggi is a matrix with 9 columns containing the annotations and information related to the signature itself.
#' @source [http://jnci.oxfordjournals.org/cgi/content/full/98/4/262/DC1](http://jnci.oxfordjournals.org/cgi/content/full/98/4/262/DC1)
#' @references Sotiriou C, Wirapati P, Loi S, Harris A, Bergh J, Smeds J, Farmer P, Praz V, Haibe-Kains B, Lallemand F, Buyse M, Piccart MJ and Delorenzi M (2006) "Gene expression profiling in breast cancer: Understanding the molecular basis of histologic grade to improve prognosis", Journal of National Cancer Institute, 98:262--272
#' @keywords data
#' @md
NULL

#' @name sig.oncotypedx
#' @docType data
#' @title Signature used to compute the OncotypeDX signature as published by Paik et al 2004
#' @description List of 21 genes included in the OncotypeDX signature. The EntrezGene.ID allows for mapping and the mapping to affy probes is already provided.
#' @usage data(sig.oncotypedx)
#' @references S. Paik, S. Shak, G. Tang, C. Kim, J. Bakker, M. Cronin, F. L. Baehner, M. G. Walker, D. Watson, T. Park, W. Hiller, E. R. Fisher, D. L. Wickerham, J. Bryant, and N. Wolmark (2004) "A Multigene Assay to Predict Recurrence of Tamoxifen-Treated, Node-Negative Breast Cancer", New England Journal of Medicine, 351(27):2817--2826.
#' @keywords data
#' @md
NULL

#' @name sig.pik3cags
#' @docType data
#' @title Gene expression Grade Index (GGI) as published in Sotiriou et al. 2006
#' @description List of 278 affymetrix hgu133a probesets representing 236 unique genes included in the PIK3CA-GS signature. The "EntrezGene.ID" column allows for mapping and "coefficient" refers to to the direction of association with PIK3CA mutation.
#' @usage data(sig.pik3cags)
#' @format sig.pik3cags is a matrix with 3 columns containing the annotations and information related to the signature itself.
#' @source [http://www.pnas.org/content/107/22/10208/suppl/DCSupplemental](http://www.pnas.org/content/107/22/10208/suppl/DCSupplemental)
#' @references Loi S, Haibe-Kains B, Majjaj S, Lallemand F, Durbecq V, Larsimont D, Gonzalez-Angulo AM, Pusztai L, Symmans FW, Bardelli A, Ellis P, Tutt AN, Gillett CE, Hennessy BT., Mills GB, Phillips WA, Piccart MJ, Speed TP, McArthur GA, Sotiriou C (2010) "PIK3CA mutations associated with gene signature of low mTORC1 signaling and better outcomes in estrogen receptor-positive breast cancer", Proceedings of the National Academy of Sciences, 107(22):10208--10213
#' @keywords data
#' @md
NULL

#' @name sig.tamr13
#' @docType data
#' @title Tamoxifen Resistance signature composed of 13 gene clusters (TAMR13) as published by Loi et al. 2008.
#' @description List of 13 clusters of genes (and annotations) and their corresponding coefficient as an additional attribute.
#' @usage data(sig.tamr13)
#' @format sig.tamr13 is a list a 13 clusters of genes with their corresponding coefficient.
#' @references Loi S, Haibe-Kains B, Desmedt C, Wirapati P, Lallemand F, Tutt AM, Gillet C, Ellis P, Ryder K, Reid JF, Daidone MG, Pierotti MA, Berns EMJJ, Jansen MPHM, Foekens JA, Delorenzi M, Bontempi G, Piccart MJ and Sotiriou C (2008) "Predicting prognosis using molecular profiling in estrogen receptor-positive breast cancer treated with tamoxifen", BMC Genomics, 9(1):239
#' @keywords data
#' @md
NULL

#' @name sigOvcAngiogenic
#' @title sigOvcAngiogenic dataset
#' @docType data
#' @source [http://jnci.oxfordjournals.org/cgi/content/full/98/4/262/DC1](http://jnci.oxfordjournals.org/cgi/content/full/98/4/262/DC1)
#' @references Bentink S, Haibe-Kains B, Risch T, Fan J-B, Hirsch MS, Holton K, Rubio R, April C, Chen J, Wickham-Garcia E, Liu J, Culhane AC, Drapkin R, Quackenbush JF, Matulonis UA (2012) "Angiogenic mRNA and microRNA Gene Expression Signature Predicts a Novel Subtype of Serous Ovarian Cancer", PloS one, 7(2):e30269
#' @keywords data
#' @md
NULL

#' @name sigOvcCrijns
#' @title sigOvcCrijns dataset
#' @docType data
#' @source [http://jnci.oxfordjournals.org/cgi/content/full/98/4/262/DC1](http://jnci.oxfordjournals.org/cgi/content/full/98/4/262/DC1)
#' @references Crijns APG, Fehrmann RSN, de Jong S, Gerbens F, Meersma G J, Klip HG, Hollema H, Hofstra RMW, te Meerman GJ, de Vries EGE, van der Zee AGJ (2009) "Survival-Related Profile, Pathways, and Transcription Factors in Ovarian Cancer" PLoS Medicine, 6(2):e1000024.
#' @keywords data
#' @md
NULL

#' @name sigOvcSpentzos
#' @title sigOcvSpentzos dataset
#' @docType data
#' @source [http://jnci.oxfordjournals.org/cgi/content/full/98/4/262/DC1](http://jnci.oxfordjournals.org/cgi/content/full/98/4/262/DC1)
#' @references Spentzos, D., Levine, D. A., Ramoni, M. F., Joseph, M., Gu, X., Boyd, J., et al. (2004). "Gene expression signature with independent prognostic significance in epithelial ovarian cancer". Journal of clinical oncology, 22(23), 4700--4710. doi:10.1200/JCO.2004.04.070
#' @keywords data
#' @md
NULL

#' @name sigOvcTCGA
#' @title sigOvcTCGA dataset
#' @docType data
#' @source [http://jnci.oxfordjournals.org/cgi/content/full/98/4/262/DC1](http://jnci.oxfordjournals.org/cgi/content/full/98/4/262/DC1)
#' @references Bell D, Berchuck A, Birrer M et al. (2011) "Integrated genomic analyses of ovarian carcinoma", Nature, 474(7353):609--615
#' @keywords data
#' @md
NULL

#' @name sigOvcYoshihara
#' @title sigOvcYoshihara dataset
#' @docType data
#' @source [http://jnci.oxfordjournals.org/cgi/content/full/98/4/262/DC1](http://jnci.oxfordjournals.org/cgi/content/full/98/4/262/DC1)
#' @references Yoshihara K, Tajima A, Yahata T, Kodama S, Fujiwara H, Suzuki M, Onishi Y, Hatae M, Sueyoshi K, Fujiwara H, Kudo, Yoshiki, Kotera K, Masuzaki H, Tashiro H, Katabuchi H, Inoue I, Tanaka K (2010) "Gene expression profile for predicting survival in advanced-stage serous ovarian cancer across two independent datasets", PloS one, 5(3):e9615.
#' @keywords data
#' @md
NULL

#' @name ssp2003
#' @aliases ssp2003.robust ssp2003.scale
#' @title SSP2003 classifier for identification of breast cancer molecular subtypes (Sorlie et al 2003)
#' @description List of parameters defining the SSP2003 classifier for identification of breast cancer molecular subtypes (Sorlie et al 2003).
#' @usage
#' data(ssp2003)
#' data(ssp2003.robust)
#' data(ssp2003.scale)
#' @docType data
#' @format List of parameters for SSP2003:
#'   - centroids: Gene expression centroids for each subtype.
#'   - centroids.map: Mapping for centroids.
#'   - method.cor: Method of correlation used to compute distance to the centroids.
#'   - method.centroids: Method used to compute the centroids.
#'   - std: Method used to compute the centroids.
#'   - mins: Minimum number of samples within each cluster allowed during the fitting of the model.
#' @source [http://www.pnas.org/content/100/14/8418](http://www.pnas.org/content/100/14/8418)
#' @references T. Sorlie and R. Tibshirani and J. Parker and T. Hastie and J. S. Marron and A. Nobel and S. Deng and H. Johnsen and R. Pesich and S. Geister and J. Demeter and C. Perou and P. E. Lonning and P. O. Brown and A. L. Borresen-Dale and D. Botstein (2003) "Repeated Observation of Breast Tumor Subtypes in Independent Gene Expression Data Sets", Proceedings of the National Academy of Sciences, 1(14):8418--8423
#' @keywords data
#' @md
NULL

#' @name ssp2006
#' @aliases ssp2006.robust ssp2006.scale
#' @title SSP2006 classifier for identification of breast cancer molecular subtypes (Hu et al 2006)
#' @description List of parameters defining the SSP2006 classifier for identification of breast cancer molecular subtypes (Hu et al 2006).
#' @usage
#' data(ssp2006)
#' data(ssp2006.robust)
#' data(ssp2006.scale)
#' @format List of parameters for SSP2006:
#'   - centroids: Gene expression centroids for each subtype.
#'   - centroids.map: Mapping for centroids.
#'   - method.cor: Method of correlation used to compute distance to the centroids.
#'   - method.centroids: Method used to compute the centroids.
#'   - std: Method of standardization for gene expressions.
#'   - mins: Minimum number of samples within each cluster allowed during the fitting of the model.
#' @details Three versions of the model are provided, each of ones differs by the gene expressions standardization method since it has an important impact on the subtype classification:
#'   - ssp2006: Use of the official centroids without scaling of the gene expressions.
#'   - ssp2006.scale: Use of the official centroids with traditional scaling of the gene expressions (see [`base::scale()`])
#'   - ssp2006.robust: Use of the official centroids with robust scaling of the gene expressions (see [`genefu::rescale()`])
#' The model `ssp2006.robust` has been shown to reach the best concordance with the traditional clinical parameters (ER IHC, HER2 IHC/FISH and histological grade). However the use of this model is recommended only when the dataset is representative of a global population of breast cancer patients (no sampling bias, the 5 subtypes should be present).
#' @docType data
#' @source [http://www.biomedcentral.com/1471-2164/7/96](http://www.biomedcentral.com/1471-2164/7/96)
#' @references Hu, Zhiyuan and Fan, Cheng and Oh, Daniel and Marron, JS and He, Xiaping and Qaqish, Bahjat and Livasy, Chad and Carey, Lisa and Reynolds, Evangeline and Dressler, Lynn and Nobel, Andrew and Parker, Joel and Ewend, Matthew and Sawyer, Lynda and Wu, Junyuan and Liu, Yudong and Nanda, Rita and Tretiakova, Maria and Orrico, Alejandra and Dreher, Donna and Palazzo, Juan and Perreard, Laurent and Nelson, Edward and Mone, Mary and Hansen, Heidi and Mullins, Michael and Quackenbush, John and Ellis, Matthew and Olopade, Olufunmilayo and Bernard, Philip and Perou, Charles (2006) "The molecular portraits of breast tumors are conserved across microarray platforms", _BMC Genomics_, *7*(96)
#' @keywords data
#' @md
NULL

#' @name vdxs
#' @aliases data.vdxs annot.vdxs demo.vdxs
#' @docType data
#' @title Gene expression, annotations and clinical data from Wang et al. 2005 and Minn et al 2007
#' @description This dataset contains (part of) the gene expression, annotations and clinical data as published in Wang et al. 2005 and Minn et al 2007.
#' @format `vdxs` is a dataset containing three matrices:
#'   - data.vdxs: Matrix containing gene expressions as measured by Affymetrix hgu133a technology (single-channel, oligonucleotides)
#'   - annot.vdxs: Matrix containing annotations of ffymetrix hgu133a microarray platform
#'   - demo.vdxs: Clinical information of the breast cancer patients whose tumors were hybridized
#' @details This dataset represent only partially the one published by Wang et al. 2005 and Minn et al 2007. Indeed only part of the patients (150) and gene expressions (966) are contained in `data.vdxs`.
#' @source
#' [http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2034](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2034)
#'
#' [http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5327](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5327)
#' @references
#' Y. Wang and J. G. Klijn and Y. Zhang and A. M. Sieuwerts and M. P. Look and F. Yang and D. Talantov and M. Timmermans and M. E. Meijer-van Gelder and J. Yu and T. Jatkoe and E. M. Berns and D. Atkins and J. A. Foekens (2005) "Gene-Expression Profiles to Predict Distant Metastasis of Lymph-Node-Negative Primary Breast Cancer", _Lancet_, *365*:671--679
#'
#' Minn, Andy J. and Gupta, Gaorav P. and Padua, David and Bos, Paula and Nguyen, Don X. and Nuyten, Dimitry and Kreike, Bas and Zhang, Yi and Wang, Yixin and Ishwaran, Hemant and Foekens, John A. and van de Vijver, Marc and Massague, Joan (2007) "Lung metastasis genes couple breast tumor size and metastatic spread", _Proceedings of the National Academy of Sciences_, *104*(16):6740--6745
#' @keywords data
#' @md
NULL

