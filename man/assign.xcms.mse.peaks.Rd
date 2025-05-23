% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assign.xcms.mse.peaks.R
\name{assign.xcms.mse.peaks}
\alias{assign.xcms.mse.peaks}
\title{assign.xcms.mse.peaks}
\usage{
assign.xcms.mse.peaks(
  ramclustObj = NULL,
  min.cor = 0.7,
  rt.sd = 2,
  min.rt.range = 0.5
)
}
\arguments{
\item{ramclustObj}{ramclustR object to annotate.}

\item{min.cor}{numeric: minimum dataset-wide correlational r-value to be considered for assignment.}

\item{rt.sd}{numeric: multiplier for the cluster standard deviation for retention times.  used for defining the retention time window to consider for fragment ion retention time assignment. retention time window is +/-  max(min.rt.range, rt.sd*ramclustObj$clrtsd)}

\item{min.rt.range}{numeric: minimum retention time range (in seconds, for xcms data) to be considered. retention time window is +/-  max(min.rt.range, rt.sd*ramclustObj$clrtsd)}
}
\value{
$M:  The inferred molecular weight of the compound giving rise to the each spectrum
}
\description{
assign.xcms.mse.peaks
}
\details{
This function is developed to assign fragment ions to ramclustR compounds after performing peak finding in the MSe (DIA MS/MS) data from centWave XCMS output.
}
\author{
Corey Broeckling
}
\concept{RAMClustR}
\concept{mass spectrometry}
\concept{metabolomics}
\concept{xcms}
