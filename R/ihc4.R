#' @title Function to compute the IHC4 prognostic score as published by 
#'   Paik et al. in 2004.
#'
#' @description
#' This function computes the prognostic score based on four measured IHC markers
#'   (ER, PGR, HER2, Ki-67), following the algorithm as published by Cuzick et al. 2011.
#'   The user has the option to either obtain just the shrinkage-adjusted IHC4 score (IHC4) 
#'   or the overall score htat also combines the clinical score (IHC4+C)
#'
#' @usage
#' ihc4(ER, PGR, HER2, Ki67,age,size,grade,node,ana,scoreWithClinical=FALSE, na.rm = FALSE)
#'
#' @param ER ER score between 0-10, calculated as (H-score/30).
#' @param PGR Progesterone Receptor score between 0-10.
#' @param HER2 Her2/neu status (0 or 1).
#' @param Ki67 Ki67 score based on percentage of positively staining malignant cells.
#' @param age patient age.
#' @param size tumor size in cm.
#' @param grade Histological grade, i.e. low (1), intermediate (2) and high (3) grade.
#' @param node Nodal status.
#' @param ana treatment with anastrozole.
#' @param scoreWithClinical TRUE to get IHC4+C score, FALSE to get just the IHC4 score.
#' @param na.rm TRUE if missing values should be removed, FALSE otherwise.
#'
#' @return
#' Shrinkage-adjusted IHC4 score or the Overall Prognostic Score based on IHC4+C 
#'   (IHC4+Clinical Score)
#'
#' @references
#' Jack Cuzick, Mitch Dowsett, Silvia Pineda, Christopher Wale, Janine Salter, Emma Quinn,
#'   Lila Zabaglo, Elizabeth Mallon, Andrew R. Green, Ian O. Ellis, Anthony Howell, Aman U. 
#'   Buzdar, and John F. Forbes (2011) "Prognostic Value of a Combined Estrogen Receptor, 
#'   Progesterone Receptor, Ki-67, and Human Epidermal Growth Factor Receptor 2 
#'   Immunohistochemical Score and Comparison with the Genomic Health Recurrence Score
#'   in Early Breast Cancer", Journal of Clinical Oncologoy, 29(32):4273–4278.
#'
#' @examples
#' # load NKI dataset
#' data(nkis)
#' # compute shrinkage-adjusted IHC4 score
#' count<-nrow(demo.nkis)
#' ihc4(ER=sample(x=1:10, size=count,replace=TRUE),PGR=sample(x=1:10, size=count,replace=TRUE),
#' HER2=sample(x=0:1,size=count,replace=TRUE),Ki67=sample(x=1:100, size=count,replace=TRUE),
#' scoreWithClinical=FALSE, na.rm=TRUE)
#'
#' # compute IHC4+C score
#' ihc4(ER=sample(x=1:10, size=count,replace=TRUE),PGR=sample(x=1:10, size=count,replace=TRUE),
#' HER2=sample(x=0:1,size=count,replace=TRUE),Ki67=sample(x=1:100, size=count,replace=TRUE),
#' age=demo.nkis[,"age"],size=demo.nkis[ ,"size"],grade=demo.nkis[ ,"grade"],node=demo.nkis[ ,"node"],
#' ana=sample(x=0:1,size=count,replace=TRUE), scoreWithClinical=TRUE, na.rm=TRUE)
#'
#' @md
#' @export
ihc4 <-
  function(ER, PGR, HER2, Ki67,age,size,grade,node,ana,scoreWithClinical=FALSE,na.rm=FALSE)
    {
    nn <- names(ER)
    if(is.null(nn)) { nn <- paste("PATIENT", 1:length(ER), sep=".") }
    names(ER) <- names(PGR) <- names(HER2) <- names(Ki67) <- nn

    cc.ix <- complete.cases(ER, PGR, HER2, Ki67)
    ER <- ER[cc.ix]
    PGR <- PGR[cc.ix]
    HER2 <- HER2[cc.ix]

    if(length(ER) != length(PGR) || length(ER) != length(HER2) || length(ER) != length(Ki67))
    {stop("ER, PGR, HER2, and Ki67 scores must have the same length!")}
    if(!all(cc.ix) & !na.rm)  { stop("NA values are present!") }

    if(!all(is.element(PGR, c("0", "1", "2","3","4","5","6","7","8","9","10")))) {
      stop("PGR scores must be between 0 and 10")
    }

    if(!all(is.element(ER, c("0", "1", "2","3","4","5","6","7","8","9","10")))) {
      stop("ER scores must be between 0 and 10")
    }

    if(!all(is.element(HER2, c("0", "1")))) {
      #if only "0" and "1" are available, map "0" -> "1" and "1" -> "3"
      stop("her2 expression must be 0 or 1!")
    }

    if(!is.numeric(Ki67)) {
      stop("Ki67 must be numeric!")
    }

    ihc4 <- 94.7 * ((0.586*HER2)-(0.100*ER)-(0.079*PGR)+(0.240*log(1 + 10 * Ki67)))

    if(scoreWithClinical==FALSE)
    {
        names(ihc4) <- nn[cc.ix]
        ihc4.score <- rep(NA, length(cc.ix))
        names(ihc4.score) <- nn
        ihc4.score[names(ihc4)] <- ihc4
        return("score"=ihc4.score)
    }

    if(scoreWithClinical==TRUE){
        size <- size[cc.ix]
        grade <- grade[cc.ix]
        node <- node[cc.ix]

        if(length(size) != length(grade) || length(grade) != length(node)) {
          stop("size, grade and lymph node stage must have the same length!")
        }
        if(!all(cc.ix) & !na.rm)  { stop("NA values are present!") }
        if(!all(is.element(grade, c("1", "2", "3")))) {
          stop("grade must be 1, 2 or 3!")
        }
        if(!all(is.element(node, c("0","1", "2", "3")))) {
          #if only "0" and "1" are available, map "0" -> "1" and "1" -> "3"
          stop("lymph node stage must be 1, 2 or 3!")
        }
        if(!is.numeric(size)) {
          stop("tumor size (cm) must be numeric!")
        }

        Ana<-ana
        N_1to3<-ifelse(node>3,0,node)
        N4<-ifelse(node==4,node,0)
        Age<-ifelse(age>65,age,0) #Age above 65 yrs
        T_1to2<-ifelse(grade<=2,grade,0)
        T_2to3<-ifelse((grade==2 | grade==3),grade,0)
        Tabove3<-ifelse(grade>3,grade,0)
        G2<-ifelse(grade==2,grade,0)
        G3<-ifelse(grade==3,grade,0)


        ClinicalScore<- 100 * ((0.417 * N_1to3) + (1.566 * N4)
                               + (0.930 * ((0.479*T_1to2)+(0.882*T_2to3)
                                  +(1.838*Tabove3)+(0.559*G2)+(0.970*G3)+(0.130*Age)-(0.149*Ana))))

        Overall<-ihc4+ClinicalScore

        names(Overall) <- nn[cc.ix]
        Overall.score <- rep(NA, length(cc.ix))
        names(Overall.score) <- nn
        Overall.score[names(Overall)] <- ihc4
        return("score"=Overall.score)
    }
  }