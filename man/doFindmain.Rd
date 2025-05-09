% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/do.findmain.R
\name{doFindmain}
\alias{doFindmain}
\title{doFindmain}
\usage{
doFindmain(
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
}
\arguments{
\item{ramclustObj}{ramclustR object to annotate.}

\item{cmpd}{integer: vector defining compound numbers to annotated.  if NULL (default), all compounds}

\item{mode}{character: "positive" or "negative"}

\item{mzabs.error}{numeric: absolute mass deviation allowd, default = 0.01}

\item{ppm.error}{numeric: ppm mass error \emph{added} to mzabs.error, default = 10}

\item{mainpkthr}{numeric.  passed to InterpretMSSpectrum::findMAIN}

\item{ads}{character: vector of allowed adducts, i.e. c("\link{M+H}+"). if NULL, default positive mode values of H+, Na+, K+, and NH4+, as monomer, dimer, and trimer, are assigned. Negative mode include "\link{M-H}-", "\link{M+Na-2H}-", "\link{M+K-2H}-", "\link{M+CH2O2-H}-" as monomer, dimer, and trimer.}

\item{nls}{character: vector of allowed neutral losses, i.e. c("\link{M+H-H2O}+").  if NULL, an extensive list derived from CAMERA's will be used.}

\item{collision.energy}{numeric: collision energy used for fragmentation. default = 25. ideally this would come from the Spectra object raw data, but this does not work well for rampled collision energy, as only the first value is pulled into the Spectra object.}
}
\value{
an updated RAMClustR object with new slots including:

$findmain - detailed output for each MS1 spectrum

$feature.table - table of all features with findmain annotations for each.

$ms1.isotopes - list of all isotopic envelopes for all clusters. when multiple isotope clusters are detected and annotated as adducts/neutral losses, all are reported in prioritization order as listed in ads and nls options.
}
\description{
Cluster annotation function: inference of 'M' - molecular weight of the compound giving rise to each spectrum - using the InterpretMSSpectrum::findMain function
}
\details{
a partially annotated ramclustR object.  base structure is that of a standard R heirarchical clustering output, with additional slots described in ramclustR documentation (?ramclustR).  New slots added after using the interpretMSSpectrum functionality include those described below.
}
\references{
Jaeger C, ... Lisec J. Compound annotation in liquid chromatography/high-resolution mass spectrometry based metabolomics: robust adduct ion determination as a prerequisite to structure prediction in electrospray ionization mass spectra. Rapid Commun Mass Spectrom. 2017 Aug 15;31(15):1261-1266. doi: 10.1002/rcm.7905. PubMed PMID: 28499062.

Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.

Broeckling CD, Ganna A, Layer M, Brown K, Sutton B, Ingelsson E, Peers G, Prenni JE. Enabling Efficient and Confident Annotation of LC-MS Metabolomics Data through MS1 Spectrum and Time Prediction. Anal Chem. 2016 Sep 20;88(18):9226-34. doi: 10.1021/acs.analchem.6b02479. Epub 2016 Sep 8. PubMed PMID: 7560453.
}
\author{
Corey Broeckling
}
\concept{RAMClustR}
\concept{clustering}
\concept{feature}
\concept{findMain}
\concept{interpretMSSpectrum}
\concept{mass spectrometry}
\concept{metabolomics}
\concept{ramclustR}
\concept{xcms}
