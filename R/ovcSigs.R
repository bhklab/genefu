`.ovcSigs` <-
function(sigs=c("bentink2012_angiogenic", "crijns2009_sig", "yoshihara2010_sig", "spentzos2011_sig", "tcga2011_sig")) {
    for(i in 1:length(sigs)) {
        sig <- read.csv(system.file(file.path("extdata", sprintf("%s.csv", sigs[i])), package="genefu"), stringsAsFactors=TRUE)
    }
}