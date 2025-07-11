#' gc.ei.create.spectra
#'
#' Cluster annotation function: inference of 'M' - molecular weight of the compound giving rise to each spectrum - using the InterpretMSSpectrum::findMain function
#'
#' @param ramclustObj ramclustR object input. 
#' @param round.digits how many digits to round mz value to.  for nist search, best to use '0'.  for other GC-EI databases, you could use a larger value, but then you should set 'mult = 1'
#' @param mult the multiplier for m/z to ensure that rounding results in the proper 'floor' mass.  NIST recommends 0.9988 to account for CH2 mass defect increase with molecular weight.
#' @details returns ramclustObj with GC-EI spectra organized into a 'Spectra' object for searching. 
#' @return $ms1.spectrum - list of all spectra for all clusters. 
#' @return $feature.table - organized table with one row per feature.    
#' @references Jaeger C, ... Lisec J. Compound annotation in liquid chromatography/high-resolution mass spectrometry based metabolomics: robust adduct ion determination as a prerequisite to structure prediction in electrospray ionization mass spectra. Rapid Commun Mass Spectrom. 2017 Aug 15;31(15):1261-1266. doi: 10.1002/rcm.7905. PubMed PMID: 28499062.
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @references Broeckling CD, Ganna A, Layer M, Brown K, Sutton B, Ingelsson E, Peers G, Prenni JE. Enabling Efficient and Confident Annotation of LC-MS Metabolomics Data through MS1 Spectrum and Time Prediction. Anal Chem. 2016 Sep 20;88(18):9226-34. doi: 10.1021/acs.analchem.6b02479. Epub 2016 Sep 8. PubMed PMID: 7560453.
#' @concept ramclustR
#' @concept RAMClustR
#' @concept metabolomics
#' @concept mass spectrometry
#' @concept clustering
#' @concept feature
#' @concept findMain
#' @concept interpretMSSpectrum
#' @concept xcms
#' @author Corey Broeckling
#' @export 


gc.ei.create.spectra <- function (
    ramclustObj = NULL, 
    round.digits = 0,
    mult = 0.9988
) 
{
  
  if(round.digits > 0 & mult != 1) {
    warning("generally when round.digits is greater than zero, mult should be set to '1'")
  }
  
  if(is.null(ramclustObj)) {
    stop("must supply ramclustObj as input.  i.e. ramclustObj = RC", '\n')
  }
  
  cmpd <- ramclustObj$cmpd
  
  spd <- DataFrame(
    name = cmpd,
    msLevel = rep(1L, length(cmpd)),
    polarity = rep(1L, length(cmpd)),
    rtime = round(ramclustObj$clrt, 1)
    )
  
  if(!is.null(ramclustObj$clri)) {
    spd$rindex = round(ramclustObj$clri, 1)
    }
  
  spd$mz <- as.list(rep(NA, length(cmpd)))
  spd$intensity <- as.list(rep(NA, length(cmpd)))
  spd$feature.index <- as.list(rep(NA, length(cmpd)))
  
  mzs <- round(mult * ramclustObj$fmz, round.digits)
  rts <- round(ramclustObj$frt, 1)
  median.intensity <- ramclustObj$msint
  feat.names <- ramclustObj$featnames

  for (cl in 1:length(cmpd)) {
    cl.feat <- which(ramclustObj$featclus == cl)
    spd$mz[[cl]] <- round(mzs[cl.feat], 1)
    mz.order <- order(spd$mz[[cl]], decreasing = FALSE)
    spd$mz[[cl]] <- spd$mz[[cl]][mz.order]
    spd$intensity[[cl]] <- median.intensity[cl.feat][mz.order]
    spd$feature.index[[cl]] <- cl.feat[mz.order]
  }
  
  ramclustObj$ms1.spectrum <- Spectra::Spectra(spd)
  # plotSpectra(ramclustObj$ms1.spectrum[1])

  # generate master findmain summary table
  feature.table <- data.frame(
    'xcms.name' = ramclustObj$featnames,
    'xcms.order' = ramclustObj$xcmsOrd,
    'feature.index' = 1:length(ramclustObj$featnames),
    'feature.mz' = round(ramclustObj$fmz, 1),
    'feature.rt' = round(ramclustObj$frt, 1), 
    'feature.med.int' = round(ramclustObj$msint),
    'rc.cluster.number' = ramclustObj$featclus
  )
  
  ramclustObj$feature.table <- feature.table
  
  return(ramclustObj)
}
