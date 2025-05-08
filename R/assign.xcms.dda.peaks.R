#' assign.xcms.dda.peaks
#'
#' 
#'
#' @param ramclustObj ramclustR object to annotate. 
#' @param xcms.feature.msmsObj object containing MS/MS spectra, generally having 
#' been created using 'xcms::featureSpectra', which should have beend derived 
#' from the same xcms object used to create the ramclustR object used here.
#' @details This function is developed to assign MS/MS spectra (generally 
#' via DDA MS/MS) from 'xcms::featureSpectra' to RAMClustR objects. Feature 
#' links are provided through XCMS feature names. 
#' this association enables export of DDA MS/MS spectra with the precursur ion 
#' isotope envelope data which can be used, for example, by Sirius.
#' @return  RAMClustR object, with a new slot called $dda.spectra, which is a 
#' list of lists.  the main list is the same length as the number of clusters, 
#' the nested list contains data for each DDA spectra which maps to an annotated 
#' mass in the ms1 spectrum for the ramclustR object.
#' @concept RAMClustR
#' @concept metabolomics
#' @concept mass spectrometry
#' @concept xcms
#' @author Corey Broeckling
#' @export 
#' 
#' 
assign.xcms.dda.peaks <- function(
    ramclustObj = NULL,
    xcms.feature.msmsObj = NULL
) {
  
  if(is.null(ramclustObj)) {
    stop("must supply ramclustObj as input.  i.e. ramclustObj = RC", '\n')
  }
  
  msms.feat.names <- xcms.feature.msmsObj$feature_id
  msms.intensity <- xcms.feature.msmsObj$totIonCurrent
  
  # add.priority <- unlist(strsplit(ramclustObj$params$findmain["ads"], ", "))
  ## change this to be ANY DDA spectrum for a compound
  dda.spectra <- as.list(rep(NA, length(ramclustObj$clrt)))
  for(i in 1:max(ramclustObj$featclus)) {
    dda.spectra[[i]] <- as.list(rep(NA, 0))
    for(j in 1:length(ramclustObj$ms1.isotopes[[i]])) {

      tar <- ramclustObj$ms1.isotopes[[i]][[j]]
      tar.feat <- tar$xcms.name
      do <- which(msms.feat.names == tar.feat)
      if(length(do) == 0) next
      if(length(do) > 1) {
        do <- do[which.max(msms.intensity[do])]
      }
      s <- data.frame(
        mz = round(unlist(xcms::mz(xcms.feature.msmsObj[do])), 5), 
        int = round(unlist(xcms::intensity(xcms.feature.msmsObj[do])), 2)
      )
      tar$spectrum.type <- "dda"
      tar$spectrum <- s
      dda.spectra[[i]][[length(dda.spectra[[i]])+1]] <- tar
    }
  }
  
  ramclustObj$dda.spectra <- dda.spectra
  return(ramclustObj)
}
