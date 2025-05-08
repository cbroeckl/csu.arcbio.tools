#' doFindmain
#'
#' Cluster annotation function: inference of 'M' - molecular weight of the compound giving rise to each spectrum - using the InterpretMSSpectrum::findMain function
#'
#' @param ramclustObj ramclustR object to annotate. 
#' @param cmpd integer: vector defining compound numbers to annotated.  if NULL (default), all compounds
#' @param mode character: "positive" or "negative"
#' @param mainpkthr numeric.  passed to InterpretMSSpectrum::findMAIN
#' @param mzabs.error numeric: absolute mass deviation allowd, default = 0.01
#' @param ppm.error numeric: ppm mass error _added_ to mzabs.error, default = 10
#' @param ads character: vector of allowed adducts, i.e. c("[M+H]+"). if NULL, default positive mode values of H+, Na+, K+, and NH4+, as monomer, dimer, and trimer, are assigned. Negative mode include "[M-H]-", "[M+Na-2H]-", "[M+K-2H]-", "[M+CH2O2-H]-" as monomer, dimer, and trimer.
#' @param nls  character: vector of allowed neutral losses, i.e. c("[M+H-H2O]+").  if NULL, an extensive list derived from CAMERA's will be used. 
#' @param collision.energy numeric: collision energy used for fragmentation. default = 25. ideally this would come from the Spectra object raw data, but this does not work well for rampled collision energy, as only the first value is pulled into the Spectra object.
#' @details a partially annotated ramclustR object.  base structure is that of a standard R heirarchical clustering output, with additional slots described in ramclustR documentation (?ramclustR).  New slots added after using the interpretMSSpectrum functionality include those described below. 
#' @return  an updated RAMClustR object with new slots including: 
#' @return $findmain - detailed output for each MS1 spectrum
#' @return $feature.table - table of all features with findmain annotations for each.  
#' @return $ms1.isotopes - list of all isotopic envelopes for all clusters. when multiple isotope clusters are detected and annotated as adducts/neutral losses, all are reported in prioritization order as listed in ads and nls options.    
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


doFindmain <- function (
    ramclustObj = NULL, 
    cmpd = NULL, 
    mode = "positive", 
    mzabs.error = 0.005, 
    ppm.error = 10,
    mainpkthr = 0.15,
    ads = NULL, 
    nls = NULL,
    collision.energy = 25
) 
{
  
  if (!requireNamespace("InterpretMSSpectrum", quietly = TRUE)) {
    stop("The use of this function requires package 'InterpretMSSpectrum'.")
  }
  
  if(is.null(ramclustObj)) {
    stop("must supply ramclustObj as input.  i.e. ramclustObj = RC", '\n')
  }
  
  
  if (is.null(ads)) {
    if (grepl("p", mode)) {
      ads <- c("[M+H]+", "[M+NH4]+", "[M+Na]+", "[M+K]+",  
               "[2M+H]+", "[2M+NH4]+", "[2M+Na]+", "[2M+K]+", 
               "[3M+H]+", "[3M+NH4]+", "[3M+Na]+", "[3M+K]+")
    }
    if (grepl("n", mode)) {
      ads <- c("[M-H]-", "[M+CH2O2-H]-", "[M+Na-2H]-", "[M+K-2H]-",  
               "[2M-H]-", "[2M+CH2O2- H]-", "[2M+Na-2H]-", "[2M+K-2H]-", 
               "[3M-H]-", "[3M+CH2O2- H]-", "[3M+Na-2H]-", "[3M+K-2H]-")
    }
    if (is.null(ads)) {
      stop("please define adducts using 'ads' or set mode to either'positive' or 'negative'")
    }
  }
  if (is.null(nls)) {
    if (grepl("p", mode)) {
      nls <- c("[M+H-H2O]+", "[M+H-NH3]+", "[M+H-HCOOH]+", 
               "[M+H-C6H12O6]+", "[M+H-C5H10O5]+", "[M+H-C12H22O11]+",
               "[M+H-COCH2]+", "[M+H-C2H3NO]+")
    }
    if (grepl("n", mode)) {
      nls <- c("[M-H-H2O]-", "[M-H-NH3]-", "[M-H-CO2]-", "[M-H-HCOOH]-", 
               "[M-H-NH3-CO2]-", "[M-H-C6H12O6]-", 
               "[M-H-C5H10O5]-", "[M-H-C12H22O11]-", "[M-H-COCH2]-" )
    }
    if (is.null(nls)) {
      stop("please define neutral losses using 'nls' or set mode to either'positive' or 'negative'")
    }
  }
  
  findmain <- as.list(rep(NA, length(ramclustObj$cmpd)))
  names(findmain) <-  ramclustObj$cmpd
  
  if(is.null(cmpd)) {
    cmpd <- (1:max(ramclustObj$featclus))
  }
  
  mzs <- ramclustObj$fmz
  rts <- ramclustObj$frt
  median.intensity <- ramclustObj$msint
  feat.names <- ramclustObj$featnames
  
  
  for (cl in cmpd) {
    cl.feat <- which(ramclustObj$featclus == cl)
    s <- data.frame(
      mz = mzs[cl.feat], 
      int = median.intensity[cl.feat],
      feature.index = cl.feat
    )
    row.names(s) <- cl.feat
    s <- s[order(s$mz),]
    
    
    ramclustObj$ms1.spectrum[[cl]] <- s
    
    out.fm <- InterpretMSSpectrum::findMAIN(
      s, 
      rules = c(ads, nls), 
      adducthyp = ads[grep("[M",ads, fixed = TRUE)], 
      ionmode = mode, 
      mzabs = mzabs.error, 
      ppm = ppm.error,
      mainpkthr = mainpkthr
    )
    
    fm.cl <- out.fm[[1]]
    fm.out <- summary(out.fm)
    fm.keep <- which(fm.out$total_score >= 0.7*max(fm.out$total_score))
    if(length(fm.keep) == 0) fm.keep <- 1
    fm.cl <- out.fm[fm.keep]
    fm.out <- fm.out[fm.keep,]
    
    out.list <- as.list(rep(NA, 2))
    names(out.list) <- c("summary", "details")
    out.list[[1]] <- fm.out
    out.list[[2]] <- fm.cl
    
    if (100 * round(cl/100, digits = 0) == cl) {
      cat(cl, "of", max(ramclustObj$featclus), "\n")
    }
    
    findmain[[cl]] <- out.list
  }
  
  ramclustObj$findmain <- findmain
  
  # generate master findmain summary table
  feature.table <- data.frame(
    'xcms.name' = ramclustObj$featnames,
    'xcms.order' = ramclustObj$xcmsOrd,
    'feature.index' = 1:length(ramclustObj$featnames),
    'feature.mz' = ramclustObj$fmz,
    'feature.rt' = ramclustObj$frt, 
    'feature.med.int' = ramclustObj$msint,
    'rc.cluster.number' = ramclustObj$featclus,
    'fm.hyp.rank' = NA,
    'fm.adduct' = NA,
    'fm.isogr' = NA,
    'fm.iso' = NA,
    'fm.charge' = NA,
    'fm.ppm' = NA,
    'fm.m' = NA,
    'fm.label' = NA,
    'fm.score' = NA,
    'fm.masses.explained' = NA,
    'fm.median.ppm' = NA
  )
  for (cl in cmpd) {
    tmp.sum <- findmain[[cl]]$summary
    tmp.det <- findmain[[cl]]$details[[1]]
    ind <- tmp.det$feature.index
    feature.table$fm.hyp.rank[ind] <- 1
    feature.table$fm.adduct[ind] <- tmp.det$adduct
    feature.table$fm.isogr[ind] <- tmp.det$isogr
    feature.table$fm.iso[ind] <- tmp.det$iso
    feature.table$fm.charge[ind] <- tmp.det$charge
    feature.table$fm.ppm[ind] <- tmp.det$ppm
    feature.table$fm.label[ind] <- tmp.det$label
    feature.table$fm.score[ind] <- tmp.sum$total_score[1]
    feature.table$fm.masses.explained[ind] <- tmp.sum$adducts_explained[1]
    feature.table$fm.median.ppm[ind] <- tmp.sum$medppm[1]
  }
  feature.table$fm.charge[which(is.na(feature.table$fm.charge))] <- 1
  feature.table$fm.m <- feature.table$feature.mz * feature.table$fm.m
  
  ramclustObj$feature.table <- feature.table
  
  ## create prioritized list of MS1 isotopes (when available)
  ## priority derived from adduct priority list order
  ## if no adduct signals in isogroup, then just assign precursor without isotopes
  
  ms1.iso <- as.list(rep(NA, length(ramclustObj$clrt)))
  for(i in 1:length(ramclustObj$clrt)) {
    ms <- feature.table[which(feature.table$rc.cluster.number == i),]
    tmp.ms1.iso <- as.list(rep(NA,0))
    
    ## loop through all adduct then neutral loss assigments, keeping only those with isogrp assigments
    all.ion.types <- c(ads, nls)
    tmp.names <- c("xcms.name", "ramclustr.cmpd", "precursor.mz", "precursor.m", "adduct.type", "retention.time", "collision.energy", "charge", "formula", "isotopes")
    
    for(ad in all.ion.types) {
      do <- which(ms$fm.adduct == ad & !is.na(ms$fm.isogr))
      if(length(do) == 0) next
      do.iso.gr <- ms$fm.isogr[do]
      tmp <- as.list(rep(NA, length(tmp.names)))
      names(tmp) <- tmp.names
      tmp$xcms.name <- ms$xcms.name[do]
      tmp$ramclustr.cmpd <- ramclustObj$cmpd[i]
      tmp$precursor.mz <- ms$feature.mz[do]
      tmp$precursor.m <- findmain[[i]]$summary[1,"neutral_mass"]
      tmp$adduct.type <- ms$fm.adduct[do]
      tmp$retention.time <- ramclustObj$clrt[i]
      tmp$collision.energy <- collision.energy
      tmp$charge <- ms$fm.charge[1]
      tmp.iso <- ms[which(ms$fm.isogr == do.iso.gr),]
      tmp.iso <- tmp.iso[order(tmp.iso$fm.iso),]
      tmp$isotopes <- tmp.iso
      tmp.ms1.iso[[length(tmp.ms1.iso)+1]] <- tmp
      rm(do)
    }
    for(ad in all.ion.types) {
      do <- which(ms$fm.adduct == ad & is.na(ms$fm.isogr))
      if(length(do) == 0) next
      tmp <- as.list(rep(NA, length(tmp.names)))
      names(tmp) <- tmp.names
      tmp$xcms.name <- ms$xcms.name[do]
      tmp$ramclustr.cmpd <- ramclustObj$cmpd[i]
      tmp$precursor.mz <- ms$feature.mz[do]
      tmp$precursor.m <- findmain[[i]]$summary[1,"neutral_mass"]
      tmp$adduct.type <- ms$fm.adduct[do]
      tmp$retention.time <- ramclustObj$clrt[i]
      tmp$collision.energy <- collision.energy
      tmp$charge <- ms$fm.charge[1]
      tmp$isotopes <- ms[do,]
      tmp.ms1.iso[[length(tmp.ms1.iso)+1]] <- tmp
      rm(do)
    }
    
    ms1.iso[[i]] <- tmp.ms1.iso 
    
  }
  
  ramclustObj$ms1.isotopes <- ms1.iso
 
  params <- c(
    mode = "positive", 
    mzabs.error = 0.005, 
    ppm.error = 10,
    mainpkthr = 0.15,
    ads = paste(ads, collapse = ', '), 
    nls = paste(nls, collapse = ', '),
    collision.energy = 25
  )
  
  ramclustObj$params$findmain <- params
  
  ramclustObj$history$do.findmain <- paste(
    " Molecular weight was inferred from in-source spectra (Broeckling 2016) using the do.findmain function, which calls the ", 
    "interpretMSSpectrum package (Jaeger 2016). ", 
    "Parameters for do.findmain were set to: ", 
    "mode = ", mode, ", mzabs.error = ", mzabs.error, ", ppm.error = ", 
    ppm.error, ", ads = ", paste(ads, collapse = " "), ", nls = ", 
    paste(nls, collapse = " "), ".", sep = "")
  cat("finished", "\n")
  return(ramclustObj)
}
