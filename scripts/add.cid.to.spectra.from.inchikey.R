
## add CID to NIST msp-derived spectral library
library(MetaboAnnotation)
library(Spectra)



## map NIST data to pubchem CID
load("R:/RSTOR-PMF/Software/db/nist.msp/nist23-hr/nist.hr.all.Rdata")

## CAS 
load("C:/Temp/20250828/cid.cas.Rdata")
nist.cas <- nist.hr.all$CAS.
nist.cas <- sapply(1:length(nist.cas), FUN = function(x) {unlist(strsplit(nist.cas[x], ";"))[1]})
## sort by cid so that the first match is the lowest cid available when multiple matches (frequent)
cid.cas.sorted <- cid.cas[order(cid.cas$cid),]
cas.match <- match(nist.cas, cid.cas.sorted$cas)
cid.from.cas <- cid.cas.sorted$cid[cas.match]
rm(cid.cas, cid.cas.sorted)
gc()

## inchikey next.
load("C:/Temp/20250828/cid.inchikey.Rdata")
nist.inchikey <- nist.hr.all$inchikey
inchikey.match <- match(nist.inchikey, cid.inchikey$inchikey)
cid.from.inchikey <- cid.inchikey$cid[inchikey.match]
length(which(is.na(cid.from.inchikey)))
check.matches <- which(!is.na(cid) & !is.na(cid.from.inchikey))
length(which(cid[check.matches] != cid.from.inchikey[check.matches]))
## quite a few disagree, but the ones that do seem to be differing CIDs with the same structure.
## at least one example of a slightly different charge state: 94287 vs 10350317
## get parent for each. 


rm(cid.inchikey)
gc()
load("C:/Temp/20250828/cid.parent.Rdata")
cas.parent.match <- match(cid, cid.parent$cid)
cas.parent <- cid.parent$parent.cid[cas.parent.match]
inchi.parent.match <- match(cid.from.inchikey, cid.parent$cid)
inchi.parent <- cid.parent$parent.cid[inchi.parent.match]
mismatch <- which(cas.parent != inchi.parent)
smp <- sample(mismatch, 100)
data.frame(cas.parent[smp], inchi.parent[smp])
## still lots of mismatches.  675014 vs 6619004 stereochemistry S/R
## 73342 vs 47641 - double bond defined vs not
## 5486187 vs 14077782 - stereochemistry defined vs not

## seems like a problem with match selecting only the first match, when there are multiple. 
## try this inchikey:  OZGNYLLQHRPOBR-DHZHZOJOSA-N and CAS: 65472-88-0
grep("OZGNYLLQHRPOBR-DHZHZOJOSA-N", cid.inchikey$inchikey)
grep("65472-88-0", cid.cas$cas)
## one match to inchikey, six matches to cas. 
## inchikey is actually better, more specific.
## inchikey also has zero missing CIDs
## skip CAS

cid <- cid.from.inchikey
cid.parent <- inchi.parent
na.rep <- which(is.na(cid.parent))
cid.parent[na.rep] <- cid[na.rep]

na.rep <- which(is.na(cas.parent))
cas.parent[na.rep] <- cid.from.cas[na.rep]

na.rep <- which(is.na(cid))
cid[na.rep] <- cid.from.cas[na.rep]
cid.parent[na.rep] <- cas.parent[na.rep]

length(which(is.na(cid)))
## 45586 missing cid
missing.cid.names <- unique(nist.hr.all$name[which(is.na(cid))])
## 1017 by name
missing.cid.names[sample(1:length(missing.cid.names), 10)]

## tried fuzzy matching, didn't like it.  Even the best match would return the wrong structure. didn't dig deeper.
## one could consider the best match which also matches formula, but it would take a good deal of evaluation to quantify acccuracy.
# install.packages('stringdist')
# library(stringdist)
# load("C:/Temp/20250828/cid.synonym.Rdata")
# syn.match <- match(missing.cid.names, cid.synonym$synonym)
# 
# matched_data <- stringdist(
#   missing.cid.names[1], cid.synonym$synonym,
#   method = "jw",  # Jaro-Winkler distance metric
# )
# which(matched_data - min(matched_data) < 0.01)
# which(matched_data == min(matched_data))
# missing.cid.names[1]
# ## not good - best match: "1-eicosanoyl-sn-glycero-3-phosphoethanolamine" for the search string "1-Eicosatrienoyl-sn-glycero-3-phosphoethanolamine"
# ## probably not worth keeping down this track.  
nist.hr.all <- Spectra::setBackend(nist.hr.all, backend = MsBackendMemory())
nist.hr.all$cid <- cid
nist.hr.all$cid.parent <- cid.parent
nist.hr.all$cas <- nist.cas

save(nist.hr.all, file = "R:/RSTOR-PMF/Software/db/nist.msp/nist23-hr/nist.hr.all.Rdata")

nist.hr.pos <- Spectra::filterPolarity(nist.hr.all, 1L)
head(sort(table(nist.hr.pos$adduct), decreasing = TRUE), 10)
save(nist.hr.pos, file = "R:/RSTOR-PMF/Software/db/nist.msp/nist23-hr/nist.hr.pos.Rdata")

nist.hr.neg <- Spectra::filterPolarity(nist.hr.all, 0L)
head(sort(table(nist.hr.neg$adduct), decreasing = TRUE), 10)
save(nist.hr.neg, file = "R:/RSTOR-PMF/Software/db/nist.msp/nist23-hr/nist.hr.neg.Rdata")


library(csu.pmf.tools)
rc.cmpd.get.pubchem(cmpd.names = missing.cid.names, 
                    get.properties = FALSE, get.bioassays = FALSE, 
                    get.vendors = FALSE, get.synonyms = FALSE, get.pathways = FALSE)



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
nist.hr.all <- Spectra::setBackend(nist.hr.all, backend = MsBackendMemory())


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
