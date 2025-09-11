library(pubchem.bio)
library(csu.arcbio.tools)
packageVersion('pubchem.bio')
local.pubchem.directory <- "C:/Temp/20250828"  ## i suggest this naming scheme, but feel free to chose your own.
load(paste0(local.pubchem.directory, "/pc.bio.Rdata"))
pc.bio <- build.element.count(pubchem.bio.object = pc.bio)
pc.bio.fiora <- pubchem.bio.to.fiora(
  pubchem.bio.object = pc.bio,
  param.table = data.frame(
    Precursor_type = c(rep("[M+H]+", 3), rep("[M-H]-", 3)), 
    CE = c(10, 20, 40, 10, 20, 40), 
    Instrument_type = rep("HCD", 6)),
  out.dir = 'C:/Temp/for.fiora/'
)




## at this point you need to run Fiora on Linux
## once complete, you can import files
## into 'Spectra' object to enable
## spectral searching


## we then import the spectra into a 'Spectra' object
mgf.dir <- "C:/Temp/20250828/fiora/"
fiora.lib <- fiora.mgf.to.spectra(mgf.dir = mgf.dir)


## try search against NIST (ONLY works from ARC-BIO computers)
library(MetaboAnnotation)
library(Spectra)
load("R:/RSTOR-PMF/Software/db/nist.msp/nist23-hr/nist.hr.pos.Rdata")
load("C:/Temp/20250828/fiora/fiora.pos.Rdata")


csp.fr <- MatchForwardReverseParam(
  MAPFUN = joinPeaks,
  tolerance = 0,
  ppm = 5,
  FUN = MsCoreUtils::ndotproduct,
  requirePrecursor = TRUE,
  requirePrecursorPeak = FALSE,
  THRESHFUN = function(x) which(x >= 0.5),
  THRESHFUN_REVERSE = NULL,
  toleranceRt = Inf,
  percentRt = 0
)


csp.dotprod <- CompareSpectraParam(
  MAPFUN = joinPeaks,
  tolerance = 0,
  ppm = 5,
  FUN = MsCoreUtils::ndotproduct,
  requirePrecursor = TRUE,
  requirePrecursorPeak = FALSE,
  THRESHFUN = function(x) which(x >= 0.5),
  toleranceRt = Inf,
  percentRt = 0
)

## not sure it will ever finish....
csp.entropy <- CompareSpectraParam(
  MAPFUN = joinPeaks,
  tolerance = 0,
  ppm = 5,
  # FUN = MsCoreUtils::nentropy,
  FUN = msentropy::msentropy_similarity,
  requirePrecursor = TRUE,
  requirePrecursorPeak = FALSE,
  THRESHFUN = function(x) which(x >= 0.5),
  toleranceRt = Inf,
  percentRt = 0
)

tests <- sample(1:length(fiora.pos), 100)
test.lib <- fiora.pos[tests]
nist.hr.pos <- Spectra::setBackend(nist.hr.pos, backend = MsBackendMemory())


a <- Sys.time()
matches.dotprod <- MetaboAnnotation::matchSpectra(
  test.lib[1:40],
  fiora.pos[1:400],
  param = csp.dotprod,
  rtColname = c("rtime", "rtime"),
  BPPARAM = BiocParallel::SnowParam(1),
  addOriginalQueryIndex = TRUE)
b <- Sys.time()
c <- Sys.time()
matches.entropy <- MetaboAnnotation::matchSpectra(
  test.lib[1:40],
  fiora.pos[1:400],
  param = csp.entropy,
  rtColname = c("rtime", "rtime"),
  BPPARAM = BiocParallel::SnowParam(1),
  addOriginalQueryIndex = TRUE)
d <- Sys.time()
e <- Sys.time()
matches.fr <- MetaboAnnotation::matchSpectra(
  test.lib[1:40],
  fiora.pos[1:400],
  param = csp.fr,
  rtColname = c("rtime", "rtime"),
  BPPARAM = BiocParallel::SnowParam(1),
  addOriginalQueryIndex = TRUE
)
f <- Sys.time()
b-a
d-c
f-e


