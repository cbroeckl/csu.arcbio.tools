% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assign.xcms.dda.peaks.R
\name{assign.xcms.dda.peaks}
\alias{assign.xcms.dda.peaks}
\title{assign.xcms.dda.peaks}
\usage{
assign.xcms.dda.peaks(ramclustObj = NULL, xcms.feature.msmsObj = NULL)
}
\arguments{
\item{ramclustObj}{ramclustR object to annotate.}

\item{xcms.feature.msmsObj}{object containing MS/MS spectra, generally having
been created using 'xcms::featureSpectra', which should have beend derived
from the same xcms object used to create the ramclustR object used here.}
}
\value{
RAMClustR object, with a new slot called $dda.spectra, which is a
list of lists.  the main list is the same length as the number of clusters,
the nested list contains data for each DDA spectra which maps to an annotated
mass in the ms1 spectrum for the ramclustR object.
}
\description{
assign.xcms.dda.peaks
}
\details{
This function is developed to assign MS/MS spectra (generally
via DDA MS/MS) from 'xcms::featureSpectra' to RAMClustR objects. Feature
links are provided through XCMS feature names.
this association enables export of DDA MS/MS spectra with the precursur ion
isotope envelope data which can be used, for example, by Sirius.
}
\author{
Corey Broeckling
}
\concept{RAMClustR}
\concept{mass spectrometry}
\concept{metabolomics}
\concept{xcms}
