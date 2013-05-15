`rorS` <-
function(data, annot, do.mapping=FALSE, mapping, verbose=FALSE) {

  ## PAM50 classification
  data(pam50)
  sbts <- intrinsic.cluster.predict(sbt.model=pam50, data=data, annot=annot, do.mapping=do.mapping, verbose=FALSE)
  mymapping <- c("mapped"=nrow(sbts$centroids.map), "total"=nrow(pam50$centroids.map))
  ## ROR-S
  rs.unscaled <- rs <- rsrisk <- rep(NA, nrow(data))
  names(rs.unscaled) <- names(rs) <- names(rsrisk) <- rownames(data)
  rst <- 0.05 * sbts$cor[ , "Basal"] + 0.12 * sbts$cor[ , "Her2"] - 0.34 * sbts$cor[ , "LumA"] + 0.23 * sbts$cor[ , "LumA"]
  rs.unscaled[names(rst)] <- rst
  ## rescale between 0 and 100
  rs <- (rs.unscaled - quantile(rs.unscaled, probs=0.025, na.rm=TRUE)) / (quantile(rs.unscaled, probs=0.975, na.rm=TRUE) - quantile(rs.unscaled, probs=0.025, na.rm=TRUE)) * 100
  rs[!is.na(rs) & rs < 0] <- 0
  rs[!is.na(rs) & rs > 100] <- 100
  rsrisk[rs < 29] <- "Low"
  rsrisk[rs >= 29 & rs < 53] <- "Intermediate"
  rsrisk[rs >= 53] <- "High"
  rsrisk <- factor(rsrisk, levels=c("Low", "Intermediate", "High"))

	return(list("score"=rs, "risk"=rsrisk, "mapping"=mymapping, "probe"=sbts$centroids.map))
}
